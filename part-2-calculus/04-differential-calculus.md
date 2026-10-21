# Chapter 4 — Differential Calculus

[← Previous](../part-1-foundations/03-trigonometry-and-periodic-functions.md) · [Table of Contents](../00-front/table-of-contents.md) · [Next →](05-integral-calculus.md)

> **Prerequisites:** [Ch. 1](../part-1-foundations/01-mathematical-language-and-algebra.md) (functions, algebra, inequalities), [Ch. 2](../part-1-foundations/02-complex-numbers.md) (complex exponentials, $j$), [Ch. 3](../part-1-foundations/03-trigonometry-and-periodic-functions.md) (trigonometric identities, radians)
> **Learning outcomes:** After this chapter you can
> - evaluate limits (algebraically, with the squeeze theorem, with L'Hôpital's rule and with series) and decide whether a function is continuous;
> - differentiate any elementary function using the definition, the product, quotient and chain rules, implicit and logarithmic differentiation;
> - compute partial and higher-order derivatives and use them to describe how a quantity responds to each of its inputs;
> - build Taylor and Maclaurin polynomials, bound their error, and use them to linearize nonlinear devices around an operating point;
> - solve rate-of-change and related-rates problems and locate optimal design values with first- and second-derivative tests;
> - quantify how component tolerances propagate to circuit performance with absolute and relative sensitivities.

Electrical engineering is built on rates. Current is the rate at which charge passes a point, $i=dq/dt$. The voltage across an inductor is proportional to the rate of change of its current, $v=L\,di/dt$; the current into a capacitor is proportional to the rate of change of its voltage, $i=C\,dv/dt$. A transistor amplifier works because a small change in input voltage produces a larger change in output current, and that ratio is a derivative (the transconductance). Designing a filter that tolerates 5 % resistors is a question about derivatives too. Differential calculus is the mathematics of *how fast* and *how sensitive*. This chapter develops it from limits, the idea underneath every derivative, through to engineering optimization and sensitivity analysis.

Throughout, angles are in **radians** unless stated otherwise. Every derivative formula for trigonometric functions in this chapter depends on that.

## 4.1 Limits and continuity

### Why limits?

Suppose a charge meter reports $q(t)$, the charge that has passed a cross-section of a wire. Over an interval $[t, t+\Delta t]$ the *average* current is

$$
\bar i = \frac{q(t+\Delta t)-q(t)}{\Delta t}.
$$

The *instantaneous* current at time $t$ is what this ratio approaches as $\Delta t$ shrinks to zero. We cannot simply set $\Delta t=0$, because that gives $0/0$. We need a precise way to say what a quantity *approaches* without ever reaching the forbidden point. That is the limit. Limits also make precise what "a capacitor voltage cannot jump" means (continuity), what a filter does "at very high frequency" (limits at infinity), and when an equation must have a root (the intermediate value theorem).

### Definition of a limit

Informally, $\lim_{x\to a} f(x)=L$ means that $f(x)$ can be made as close to $L$ as we like by taking $x$ close enough to $a$ (but not equal to $a$). The value $f(a)$ itself, or whether it exists, is irrelevant.

> [!definition] Limit ($\varepsilon$–$\delta$)
> Let $f$ be defined on an open interval containing $a$, except possibly at $a$. We write
> $$
> \lim_{x\to a} f(x) = L
> $$
> if for every $\varepsilon>0$ there exists $\delta>0$ such that
> $$
> 0<|x-a|<\delta \quad\Longrightarrow\quad |f(x)-L|<\varepsilon .
> $$

Read it as a game. An adversary names a tolerance $\varepsilon$ on the output. You must reply with a tolerance $\delta$ on the input that guarantees the output stays within $\varepsilon$ of $L$. If you can always win, the limit is $L$. This is the same logic as a specification: "the output must be within $\pm 1\,\text{mV}$ (that is $\varepsilon$), so how tightly must the input be held (that is $\delta$)?"

> [!example] Worked Example 4.1 — An $\varepsilon$–$\delta$ proof
> Prove that $\lim_{x\to 2}(3x+1)=7$.
>
> **Step 1 (scratch work).** We need $|(3x+1)-7|<\varepsilon$, i.e. $|3x-6|=3|x-2|<\varepsilon$, i.e. $|x-2|<\varepsilon/3$.
>
> **Step 2 (choose $\delta$).** Take $\delta=\varepsilon/3$.
>
> **Step 3 (verify).** If $0<|x-2|<\delta$ then $|(3x+1)-7| = 3|x-2| < 3\delta=\varepsilon$. $\blacksquare$
>
> **Check.** With $\varepsilon=0.03$, $\delta=0.01$: $x=2.009$ gives $3x+1=7.027$, which is within $0.03$ of $7$. ✓

Nobody evaluates practical limits this way. The definition matters because every limit law below is proved from it, and because it pins down exactly what "approaches" means when intuition fails.

### One-sided limits

The **right-hand limit** $\lim_{x\to a^+}f(x)$ only lets $x$ approach from above ($x>a$); the **left-hand limit** $\lim_{x\to a^-}f(x)$ only from below. The two-sided limit exists if and only if both one-sided limits exist and are equal.

The unit step $u(t)$ (equal to $0$ for $t<0$ and $1$ for $t>0$) has $\lim_{t\to 0^-}u(t)=0$ and $\lim_{t\to 0^+}u(t)=1$. These differ, so $\lim_{t\to0}u(t)$ does not exist. Circuit analysts use one-sided limits constantly: $v_C(0^-)$ is the capacitor voltage just before a switch operates and $v_C(0^+)$ the value just after. The principle that a capacitor voltage is continuous, $v_C(0^+)=v_C(0^-)$, is a statement about one-sided limits.

### Limit laws

> [!theorem] Limit laws
> If $\lim_{x\to a}f(x)=L$ and $\lim_{x\to a}g(x)=M$ (both finite), then
> $$
> \begin{aligned}
> \lim_{x\to a}[f(x)\pm g(x)] &= L\pm M, & \lim_{x\to a}[c\,f(x)] &= cL,\\
> \lim_{x\to a}[f(x)\,g(x)] &= LM, & \lim_{x\to a}\frac{f(x)}{g(x)} &= \frac{L}{M}\quad (M\neq0),
> \end{aligned}
> $$
> and $\lim_{x\to a}[f(x)]^{n}=L^{n}$ for any positive integer $n$.

*Proof sketch of the sum law.* Given $\varepsilon$, choose $\delta_1$ so that $|f-L|<\varepsilon/2$ and $\delta_2$ so that $|g-M|<\varepsilon/2$. With $\delta=\min(\delta_1,\delta_2)$, the triangle inequality gives $|(f+g)-(L+M)|\le|f-L|+|g-M|<\varepsilon$. The product law uses $fg-LM = f(g-M)+M(f-L)$ in the same way.

Because polynomials are built from sums and products of $x$ and constants, the laws give $\lim_{x\to a}p(x)=p(a)$ for every polynomial, and $\lim_{x\to a}p(x)/q(x)=p(a)/q(a)$ for rational functions whenever $q(a)\neq0$. The interesting cases are the ones where direct substitution gives an **indeterminate form** such as $0/0$ or $\infty/\infty$.

> [!example] Worked Example 4.2 — Removing a $0/0$ by factoring
> Evaluate $\displaystyle \lim_{x\to 2}\frac{x^2-4}{x-2}$ and $\displaystyle\lim_{h\to 0}\frac{\sqrt{9+h}-3}{h}$.
>
> **First limit.** Substitution gives $0/0$. Factor: $\dfrac{(x-2)(x+2)}{x-2}=x+2$ for $x\neq 2$. Since the limit ignores the point $x=2$ itself, the limit is $\lim_{x\to2}(x+2)=4$.
>
> **Second limit.** Multiply by the conjugate:
> $$
> \frac{\sqrt{9+h}-3}{h}\cdot\frac{\sqrt{9+h}+3}{\sqrt{9+h}+3}=\frac{(9+h)-9}{h(\sqrt{9+h}+3)}=\frac{1}{\sqrt{9+h}+3}\;\longrightarrow\;\frac16 .
> $$
>
> **Numerical check.** $h=10^{-4}$: $(\sqrt{9.0001}-3)/10^{-4}=0.166666\ldots$ ✓. (This second limit is the derivative of $\sqrt x$ at $x=9$; see §4.2.)

### The squeeze theorem and the fundamental trigonometric limit

> [!theorem] Squeeze (sandwich) theorem
> If $g(x)\le f(x)\le h(x)$ for all $x$ near $a$ (except possibly at $a$) and $\lim_{x\to a}g(x)=\lim_{x\to a}h(x)=L$, then $\lim_{x\to a}f(x)=L$.

The most important application in this book is

$$
\lim_{x\to 0}\frac{\sin x}{x}=1 .
$$

*Proof.* For $0<x<\pi/2$, draw the unit circle with an angle $x$. The triangle with vertices $(0,0)$, $(1,0)$, $(\cos x,\sin x)$ has area $\tfrac12\sin x$; the circular sector has area $\tfrac12 x$; the larger right triangle with legs $1$ and $\tan x$ has area $\tfrac12\tan x$. These regions are nested, so

$$
\tfrac12\sin x\le\tfrac12 x\le\tfrac12\tan x
\quad\Longrightarrow\quad
\cos x\le\frac{\sin x}{x}\le 1 .
$$

As $x\to0^+$, $\cos x\to1$, so the squeeze theorem gives the right-hand limit $1$. Because $\sin x/x$ is an even function, the left-hand limit is also $1$. $\blacksquare$

A consequence used repeatedly below:

$$
\lim_{x\to0}\frac{1-\cos x}{x}=\lim_{x\to0}\frac{\sin^2 x}{x(1+\cos x)}=\lim_{x\to0}\frac{\sin x}{x}\cdot\frac{\sin x}{1+\cos x}=1\cdot 0=0,
$$

and in the same way $\lim_{x\to0}(1-\cos x)/x^2=\tfrac12$.

```desmos-graph
left=-15; right=15; bottom=-0.4; top=1.2
---
y=\frac{\sin(x)}{x}|blue
y=\frac{1}{x}|red|dashed|0.5<x<15
y=-\frac{1}{x}|red|dashed|0.5<x<15
(0,1)|label:limit 1
```
*The sinc-type function $\sin x/x$ is undefined at $x=0$ but tends to $1$ there; for large $x$ it is squeezed between $\pm 1/x$ and so tends to $0$.*

> [!tip] Engineering Note — the sinc function
> The function $\operatorname{sinc}(x)=\sin(\pi x)/(\pi x)$ (with $\operatorname{sinc}(0)=1$ by the limit above) is the Fourier transform of a rectangular pulse ([Ch. 16](../part-5-transforms/16-fourier-transforms.md)) and the ideal interpolation kernel in sampling theory. Defining $\operatorname{sinc}(0)=1$ is not a convention pulled from the air; it is the unique value that makes the function continuous.

> [!example] Worked Example 4.3 — Trigonometric limits
> Evaluate (a) $\displaystyle\lim_{x\to0}\frac{\sin 3x}{x}$, (b) $\displaystyle\lim_{x\to0}\frac{\tan 5x}{\sin 2x}$, (c) $\displaystyle\lim_{x\to0}\frac{1-\cos 4x}{x^2}$.
>
> (a) Write $\dfrac{\sin 3x}{x}=3\cdot\dfrac{\sin 3x}{3x}$. With $u=3x\to0$, the limit is $3\cdot1=3$.
>
> (b) $\dfrac{\tan 5x}{\sin 2x}=\dfrac{\sin 5x}{5x}\cdot\dfrac{2x}{\sin 2x}\cdot\dfrac{1}{\cos 5x}\cdot\dfrac{5}{2}\to 1\cdot1\cdot1\cdot\tfrac52=\tfrac52$.
>
> (c) With $u=4x$, $\dfrac{1-\cos u}{x^2}=16\,\dfrac{1-\cos u}{u^2}\to16\cdot\tfrac12=8$.
>
> **Numerical check** at $x=10^{-3}$: (a) $2.9999955$, (b) $2.50002$, (c) $7.99999$. ✓

### Limits at infinity and infinite limits

$\lim_{x\to\infty}f(x)=L$ means $f(x)$ can be made arbitrarily close to $L$ by taking $x$ large enough; the line $y=L$ is then a **horizontal asymptote**. $\lim_{x\to a}f(x)=\infty$ means $f(x)$ exceeds any bound as $x\to a$; the line $x=a$ is then a **vertical asymptote**. "$=\infty$" here is shorthand for a particular way of *not* having a finite limit.

For rational functions, divide numerator and denominator by the highest power of $x$ in the denominator. The useful facts are

$$
\lim_{x\to\infty}\frac{1}{x^p}=0\ (p>0),\qquad
\lim_{x\to\infty}x^n e^{-x}=0\ (\text{any }n),\qquad
\lim_{x\to\infty}\frac{\ln x}{x^p}=0\ (p>0).
$$

In words: exponentials beat powers, and powers beat logarithms. We prove the last two with L'Hôpital's rule in §4.2.

The number $e$ is itself a limit:

$$
e=\lim_{n\to\infty}\left(1+\frac1n\right)^{n}=2.718281828\ldots,
\qquad\text{and more generally}\qquad
\lim_{n\to\infty}\left(1+\frac{x}{n}\right)^{n}=e^{x}.
$$

> [!example] Worked Example 4.4 — Asymptotic behaviour of a transfer function
> A first-order low-pass filter has magnitude response
> $$
> |H(\omega)|=\frac{1}{\sqrt{1+(\omega/\omega_c)^2}} .
> $$
> Find $\lim_{\omega\to0}|H|$, $\lim_{\omega\to\infty}|H|$, and $\lim_{\omega\to\infty}\omega\,|H(\omega)|$.
>
> **Low frequency.** Direct substitution: $|H(0)|=1$. The function is continuous, so the limit is $1$ (the filter passes DC).
>
> **High frequency.** The denominator grows without bound, so $|H|\to0$.
>
> **The product.** Divide inside the root by $\omega^2$:
> $$
> \omega|H|=\frac{\omega}{\sqrt{1+\omega^2/\omega_c^2}}=\frac{1}{\sqrt{1/\omega^2+1/\omega_c^2}}\;\longrightarrow\;\omega_c .
> $$
> So for $\omega\gg\omega_c$, $|H|\approx\omega_c/\omega$: the response falls in proportion to $1/\omega$, which is $-20$ dB per decade on a Bode plot.
>
> **Check.** $\omega=1000\,\omega_c$: $\omega|H|/\omega_c = 1000/\sqrt{1+10^6}=0.9999995$. ✓

### Continuity

> [!definition] Continuity
> $f$ is **continuous at $a$** if (i) $f(a)$ is defined, (ii) $\lim_{x\to a}f(x)$ exists, and (iii) $\lim_{x\to a}f(x)=f(a)$.
> $f$ is continuous on an interval if it is continuous at every point of it (using one-sided limits at closed endpoints).

Sums, products, quotients (where the denominator is nonzero) and compositions of continuous functions are continuous. Polynomials, rational functions (on their domains), $e^x$, $\ln x$ ($x>0$), $\sin x$, $\cos x$, and $\sqrt x$ ($x\ge0$) are all continuous. So nearly every formula you meet is continuous wherever it is defined, and the action is at the points where it is not.

**Types of discontinuity:**

| Type | What fails | Example |
|---|---|---|
| Removable | limit exists but $\ne f(a)$ or $f(a)$ undefined | $\sin x/x$ at $0$ |
| Jump | one-sided limits exist but differ | $u(t)$ at $0$; a switched source |
| Infinite | function unbounded near $a$ | $1/(x-1)$ at $x=1$; an undamped resonance $1/(\omega_0^2-\omega^2)$ |
| Oscillatory | no limit, no blow-up | $\sin(1/x)$ at $0$ |

```desmos-graph
left=-3; right=4; bottom=-2; top=4
---
y=x+1|blue|x<1
y=4-x|blue|x>1
(1,2)|label:left limit 2
(1,3)|red|label:right limit 3
```
*A jump discontinuity at $x=1$: the left-hand limit is $2$ and the right-hand limit is $3$, so no choice of $f(1)$ can make the function continuous there.*

> [!example] Worked Example 4.5 — Choosing a constant for continuity
> For what value of $k$ is
> $$
> f(x)=\begin{cases}\dfrac{e^{2x}-1}{x}, & x\neq0\\[2mm] k, & x=0\end{cases}
> $$
> continuous at $0$?
>
> We need $k=\lim_{x\to0}(e^{2x}-1)/x$. With $u=2x$, this is $2\lim_{u\to0}(e^u-1)/u$. The limit $\lim_{u\to0}(e^u-1)/u=1$ is the defining property of $e$ (it says that the slope of $e^u$ at $u=0$ is $1$; see §4.2). Hence $k=2$.
>
> **Check.** $x=10^{-5}$: $(e^{2\times10^{-5}}-1)/10^{-5}=2.00002$. ✓

### Two theorems about continuous functions

> [!theorem] Intermediate value theorem (IVT)
> If $f$ is continuous on $[a,b]$ and $N$ lies strictly between $f(a)$ and $f(b)$, then there is at least one $c\in(a,b)$ with $f(c)=N$. In particular, if $f(a)$ and $f(b)$ have opposite signs, $f$ has a root in $(a,b)$.

> [!theorem] Extreme value theorem (EVT)
> If $f$ is continuous on a closed, bounded interval $[a,b]$, then $f$ attains an absolute maximum and an absolute minimum on $[a,b]$.

The IVT is the basis of the **bisection method** for solving equations ([Ch. 22](../part-7-numerical-computational/22-numerical-methods.md)): if $f(a)f(b)<0$, evaluate $f$ at the midpoint and keep whichever half still has a sign change. Each step halves the interval, so $n$ steps shrink an interval of length $L$ to $L/2^n$. The EVT guarantees that the optimization problems of §4.10 have answers when the design variable lives on a closed interval.

> [!example] Worked Example 4.6 — Bracketing a root
> Show that $x^3-x-1=0$ has a root in $(1,2)$ and locate it to within $0.1$.
>
> $f(x)=x^3-x-1$ is a polynomial, hence continuous. $f(1)=-1<0$ and $f(2)=5>0$, so by the IVT there is a root in $(1,2)$.
>
> Bisection: $f(1.5)=0.875>0\Rightarrow(1,1.5)$; $f(1.25)=-0.297<0\Rightarrow(1.25,1.5)$; $f(1.375)=0.224>0\Rightarrow(1.25,1.375)$; $f(1.3125)=-0.0515<0\Rightarrow(1.3125,1.375)$. The interval now has length $0.0625<0.1$; the root is $1.34\pm0.03$.
>
> **Check.** The root to 9 digits is $1.324717957$, inside the final bracket. ✓

> [!warning] Common Pitfall — the IVT needs continuity on the whole interval
> $f(x)=1/x$ has $f(-1)=-1$ and $f(1)=1$, yet there is no root in $(-1,1)$. The IVT does not apply because $f$ is discontinuous at $0$. A numerical root finder that brackets a sign change across a pole will happily "converge" to the pole. Always check for singularities (for example, poles of a transfer function) inside a bracket.

> [!tip] Engineering Note — continuity of state variables
> In a circuit with finite voltages and currents, capacitor voltages and inductor currents are continuous functions of time. The reason is that $i=C\,dv/dt$ is finite, so $v$ cannot jump (a jump would need infinite current); dually for $v=L\,di/dt$. This continuity, $v_C(0^+)=v_C(0^-)$ and $i_L(0^+)=i_L(0^-)$, supplies the initial conditions for every switched-circuit problem in [Ch. 11](../part-4-differential-equations/11-ordinary-differential-equations.md) and Ch. 12. Currents through capacitors and voltages across inductors are *not* state variables and may jump.

## 4.2 Derivatives and differentiation rules

### The derivative as a limit

> [!definition] Derivative
> The **derivative** of $f$ at $x$ is
> $$
> f'(x)=\lim_{h\to0}\frac{f(x+h)-f(x)}{h},
> $$
> provided the limit exists; $f$ is then **differentiable** at $x$. Equivalent notations: $f'(x)$, $\dfrac{df}{dx}$, $\dfrac{d}{dx}f(x)$, $Df(x)$, and $\dot f$ when the variable is time.

The **difference quotient** $\frac{f(x+h)-f(x)}{h}$ is the slope of the secant line through $(x,f(x))$ and $(x+h,f(x+h))$. As $h\to0$ the secant rotates into the **tangent line**, whose equation at $x=a$ is

$$
y=f(a)+f'(a)(x-a).
$$

Physically the derivative is an **instantaneous rate of change**, with units of (units of $f$)/(units of $x$). If $q$ is in coulombs and $t$ in seconds, $dq/dt$ is in C/s = A.

```desmos-graph
left=-0.5; right=3.5; bottom=-0.5; top=5
---
y=0.5x^2+0.5|blue
y=1.5x-0.625|red
y=2.25x-1.75|grey|dashed
(1.5,1.625)|label:tangent point
(2.5,3.625)|label:secant point
```
*The curve $y=\tfrac12x^2+\tfrac12$, its tangent at $x=1.5$ (slope $1.5$), and a secant through $x=1.5$ and $x=2.5$ (slope $2$). As the second point slides toward the first, the secant slope tends to the tangent slope.*

> [!example] Worked Example 4.7 — Derivatives from first principles
> Find $f'(x)$ from the definition for (a) $f(x)=x^2$, (b) $f(x)=1/x$, (c) $f(x)=\sqrt x$.
>
> (a) $\dfrac{(x+h)^2-x^2}{h}=\dfrac{2xh+h^2}{h}=2x+h\to 2x$.
>
> (b) $\dfrac{1}{h}\left(\dfrac1{x+h}-\dfrac1x\right)=\dfrac{x-(x+h)}{h\,x(x+h)}=\dfrac{-1}{x(x+h)}\to-\dfrac1{x^2}$.
>
> (c) Rationalize as in Example 4.2: $\dfrac{\sqrt{x+h}-\sqrt x}{h}=\dfrac{1}{\sqrt{x+h}+\sqrt x}\to\dfrac1{2\sqrt x}$ (for $x>0$). At $x=9$ this is $1/6$, the limit computed in Example 4.2.
>
> All three agree with the power rule $\frac{d}{dx}x^n=nx^{n-1}$ for $n=2,-1,\tfrac12$.

### Differentiability and continuity

> [!theorem] Differentiable implies continuous
> If $f$ is differentiable at $a$, then $f$ is continuous at $a$.

*Proof.* $f(a+h)-f(a)=\dfrac{f(a+h)-f(a)}{h}\cdot h\to f'(a)\cdot0=0$. $\blacksquare$

The converse is false. $f(x)=|x|$ is continuous at $0$, but its difference quotient $|h|/h$ is $+1$ for $h>0$ and $-1$ for $h<0$, so the one-sided limits disagree and $f'(0)$ does not exist. Graphically, a **corner** (or a vertical tangent, as for $x^{1/3}$ at $0$, or a cusp) blocks differentiability. An ideal diode's $i$–$v$ characteristic and the output of a full-wave rectifier $|\sin\omega t|$ both have corners. At a corner the left and right derivatives, defined with one-sided limits, exist but differ.

### Derivatives of the basic functions

**Powers.** For a positive integer $n$, the binomial theorem gives $(x+h)^n=x^n+nx^{n-1}h+\binom n2x^{n-2}h^2+\cdots$, so

$$
\frac{(x+h)^n-x^n}{h}=nx^{n-1}+\binom n2x^{n-2}h+\cdots+h^{n-1}\to nx^{n-1}.
$$

The rule $\frac d{dx}x^r=rx^{r-1}$ holds for every real $r$ (for $x>0$); we prove the general case with logarithmic differentiation in §4.3.

**Sine and cosine.** Using $\sin(x+h)=\sin x\cos h+\cos x\sin h$,

$$
\frac{\sin(x+h)-\sin x}{h}=\sin x\,\frac{\cos h-1}{h}+\cos x\,\frac{\sin h}{h}\;\longrightarrow\;\sin x\cdot0+\cos x\cdot1=\cos x.
$$

A similar computation gives $\frac d{dx}\cos x=-\sin x$. Both rely on the two trigonometric limits from §4.1, which in turn rely on measuring $x$ in radians.

**Exponential.** $\dfrac{e^{x+h}-e^x}{h}=e^x\,\dfrac{e^h-1}{h}\to e^x$, because $e$ is precisely the base for which $\lim_{h\to0}(e^h-1)/h=1$. The exponential is the function equal to its own derivative, which is why it solves every first-order linear circuit.

**Logarithm.** From $\ln(x+h)-\ln x=\ln(1+h/x)$ and the substitution $n=x/h$,

$$
\frac{\ln(1+h/x)}{h}=\frac1x\ln\left(1+\frac{1}{n}\right)^{n}\to\frac1x\ln e=\frac1x .
$$

### The differentiation rules

> [!theorem] Rules of differentiation
> If $f$ and $g$ are differentiable and $c$ is a constant:
> $$
> \begin{aligned}
> &\text{Linearity:} && (cf\pm g)'=cf'\pm g'\\
> &\text{Product rule:} && (fg)'=f'g+fg'\\
> &\text{Quotient rule:} && \left(\frac fg\right)'=\frac{f'g-fg'}{g^2}\quad(g\neq0)\\
> &\text{Power rule:} && (x^r)'=rx^{r-1}
> \end{aligned}
> $$

*Proof of the product rule.* Add and subtract $f(x+h)g(x)$:

$$
\frac{f(x+h)g(x+h)-f(x)g(x)}{h}=f(x+h)\frac{g(x+h)-g(x)}{h}+g(x)\frac{f(x+h)-f(x)}{h}.
$$

As $h\to0$, $f(x+h)\to f(x)$ (differentiable implies continuous), and the two quotients tend to $g'$ and $f'$. $\blacksquare$

*Quotient rule.* Write $f=(f/g)\,g$ and apply the product rule: $f'=(f/g)'g+(f/g)g'$. Solve for $(f/g)'$.

A useful mnemonic for the product rule is in differentials: $d(fg)=g\,df+f\,dg$. If $f$ and $g$ each change by a small amount, the product changes by "each one's change times the other". For a capacitor with charge $q=Cv$ whose capacitance *and* voltage vary, this gives $i=\dot q=\dot C v+C\dot v$ (see §4.9).

### Table of standard derivatives

| $f(x)$ | $f'(x)$ | $f(x)$ | $f'(x)$ |
|---|---|---|---|
| $x^r$ | $rx^{r-1}$ | $e^{ax}$ | $ae^{ax}$ |
| $a^x$ | $a^x\ln a$ | $\ln\lvert x\rvert$ | $1/x$ |
| $\log_a x$ | $1/(x\ln a)$ | $\sin x$ | $\cos x$ |
| $\cos x$ | $-\sin x$ | $\tan x$ | $\sec^2x=1+\tan^2x$ |
| $\cot x$ | $-\csc^2 x$ | $\sec x$ | $\sec x\tan x$ |
| $\arcsin x$ | $1/\sqrt{1-x^2}$ | $\arccos x$ | $-1/\sqrt{1-x^2}$ |
| $\arctan x$ | $1/(1+x^2)$ | $\sinh x$ | $\cosh x$ |
| $\cosh x$ | $\sinh x$ | $\tanh x$ | $\operatorname{sech}^2 x=1-\tanh^2 x$ |

The inverse-trigonometric entries follow from the inverse-function rule in §4.3. The hyperbolic ones follow from $\sinh x=(e^x-e^{-x})/2$, $\cosh x=(e^x+e^{-x})/2$.

**Complex exponentials.** For a real variable $t$ and a complex constant $s$, $\frac{d}{dt}e^{st}=se^{st}$. The proof is identical to the real case (or follows by differentiating $e^{\sigma t}\cos\omega t$ and $e^{\sigma t}\sin\omega t$ separately). With $s=j\omega$:

$$
\frac{d}{dt}e^{j\omega t}=j\omega\,e^{j\omega t}.
$$

This single line is why **differentiation becomes multiplication by $j\omega$** in phasor analysis ([Ch. 2](../part-1-foundations/02-complex-numbers.md)) and multiplication by $s$ in the Laplace domain ([Ch. 14](../part-5-transforms/14-laplace-transforms.md)). An inductor's $v=L\,di/dt$ becomes $V=j\omega L\,I$, giving the impedance $j\omega L$.

> [!example] Worked Example 4.8 — Product and quotient rules: power into a load
> A source of EMF $V$ and internal resistance $R_s$ drives a load $R_L$. The load power is
> $$
> P(R_L)=\frac{V^2R_L}{(R_s+R_L)^2}.
> $$
> Find $dP/dR_L$.
>
> **Quotient rule** with numerator $N=V^2R_L$ ($N'=V^2$) and denominator $D=(R_s+R_L)^2$ ($D'=2(R_s+R_L)$, by the chain rule of §4.3 or by expanding):
> $$
> \frac{dP}{dR_L}=\frac{V^2(R_s+R_L)^2-V^2R_L\cdot2(R_s+R_L)}{(R_s+R_L)^4}
> =\frac{V^2\left[(R_s+R_L)-2R_L\right]}{(R_s+R_L)^3}
> =\frac{V^2(R_s-R_L)}{(R_s+R_L)^3}.
> $$
>
> **Interpretation.** The derivative is positive for $R_L<R_s$ and negative for $R_L>R_s$: power rises, peaks at $R_L=R_s$, then falls. This is the maximum power transfer theorem, revisited in §4.10.
>
> **Units check.** $V^2\,\Omega/\Omega^3=V^2/\Omega^2=\text{W}/\Omega$. ✓
>
> **Numerical check.** $V=10\,$V, $R_s=50\,\Omega$, $R_L=25\,\Omega$: formula gives $100\cdot25/75^3=5.926\times10^{-3}$ W/Ω; central difference with $h=10^{-3}\,\Omega$ gives $5.926\times10^{-3}$. ✓

> [!example] Worked Example 4.9 — Inductor voltage
> The current in a $L=10\,$mH inductor is $i(t)=I_0\,(t/\tau)\,e^{-t/\tau}$ with $I_0=2\,$A and $\tau=1\,$ms. Find $v(t)=L\,di/dt$ and the time at which $v=0$.
>
> **Product rule.** With $f=t/\tau$ and $g=e^{-t/\tau}$, $f'=1/\tau$ and $g'=-e^{-t/\tau}/\tau$:
> $$
> \frac{di}{dt}=\frac{I_0}{\tau}e^{-t/\tau}-\frac{I_0 t}{\tau^2}e^{-t/\tau}=\frac{I_0}{\tau}\left(1-\frac t\tau\right)e^{-t/\tau}.
> $$
> So $v(t)=\dfrac{LI_0}{\tau}\left(1-\dfrac t\tau\right)e^{-t/\tau}=20\,(1-t/\tau)\,e^{-t/\tau}$ V.
>
> **Zero.** $v=0$ at $t=\tau=1\,$ms, which is exactly when the current reaches its peak $I_0/e=0.736\,$A. At $t=0$, $v=20\,$V.
>
> **Units.** $\text{H}\cdot\text{A}/\text{s}=\text{V}$. ✓

### The mean value theorem

> [!theorem] Rolle's theorem and the mean value theorem (MVT)
> Let $f$ be continuous on $[a,b]$ and differentiable on $(a,b)$.
> - (Rolle) If $f(a)=f(b)$, there is $c\in(a,b)$ with $f'(c)=0$.
> - (MVT) There is $c\in(a,b)$ with
> $$
> f'(c)=\frac{f(b)-f(a)}{b-a}.
> $$

*Proof sketch.* Rolle: by the EVT $f$ has a max and a min on $[a,b]$; if both are at the endpoints, $f$ is constant and $f'\equiv0$; otherwise an interior extremum has $f'=0$ (Fermat's theorem, §4.10). MVT: apply Rolle to $g(x)=f(x)-\frac{f(b)-f(a)}{b-a}(x-a)$, which has $g(a)=g(b)=f(a)$.

In words: somewhere, the instantaneous rate equals the average rate. If a capacitor voltage rises from $0$ to $5\,$V in $2\,\mu$s, then at some instant $dv/dt=2.5\,\text{V}/\mu\text{s}$ exactly, and with $C=100\,$pF the current at that instant is $i=C\,dv/dt=0.25\,$mA, whatever the shape of the waveform.

Three consequences used throughout the book:

1. If $f'(x)=0$ on an interval, $f$ is constant there.
2. If $f'=g'$ on an interval, $f-g$ is constant (the basis of "$+C$" in [Ch. 5](05-integral-calculus.md)).
3. If $f'>0$ on an interval, $f$ is strictly increasing there; if $f'<0$, strictly decreasing.

The MVT also gives error bounds: $|f(b)-f(a)|\le M|b-a|$ where $M=\max|f'|$. A sensor with transfer function slope at most $M$ amplifies an input uncertainty $\Delta x$ to at most $M\Delta x$.

### L'Hôpital's rule

> [!theorem] L'Hôpital's rule
> Suppose $f$ and $g$ are differentiable near $a$ (except possibly at $a$), $g'(x)\neq0$ near $a$, and $\lim f=\lim g=0$ or $\lim f=\pm\lim g=\pm\infty$. If $\lim_{x\to a}f'(x)/g'(x)=L$ (finite or infinite), then
> $$
> \lim_{x\to a}\frac{f(x)}{g(x)}=L .
> $$
> The same holds for one-sided limits and for $a=\pm\infty$.

*Why it works (for $0/0$ with $f'$, $g'$ continuous and $g'(a)\neq 0$).* Since $f(a)=g(a)=0$,

$$
\frac{f(x)}{g(x)}=\frac{f(x)-f(a)}{g(x)-g(a)}=\frac{[f(x)-f(a)]/(x-a)}{[g(x)-g(a)]/(x-a)}\to\frac{f'(a)}{g'(a)} .
$$

The general case uses Cauchy's version of the MVT.

Other indeterminate forms ($0\cdot\infty$, $\infty-\infty$, $1^\infty$, $0^0$, $\infty^0$) are first rewritten as $0/0$ or $\infty/\infty$, usually by moving a factor to the denominator or by taking logarithms.

> [!example] Worked Example 4.10 — L'Hôpital's rule
> (a) $\displaystyle\lim_{x\to0}\frac{e^x-1-x}{x^2}$. Form $0/0$. Differentiate top and bottom: $\dfrac{e^x-1}{2x}$, still $0/0$. Again: $\dfrac{e^x}{2}\to\dfrac12$.
>
> (b) $\displaystyle\lim_{x\to\infty}x^3e^{-x}=\lim\frac{x^3}{e^x}$. Form $\infty/\infty$. Three applications: $\dfrac{3x^2}{e^x}\to\dfrac{6x}{e^x}\to\dfrac{6}{e^x}\to0$. The same argument works for any power $x^n$: exponentials dominate.
>
> (c) $\displaystyle\lim_{x\to0^+}x\ln x$. Form $0\cdot(-\infty)$. Rewrite as $\dfrac{\ln x}{1/x}$ ($-\infty/\infty$): $\dfrac{1/x}{-1/x^2}=-x\to0$.
>
> (d) $\displaystyle\lim_{x\to\infty}\left(1+\frac{a}{x}\right)^x$. Form $1^\infty$. Take logs: $x\ln(1+a/x)=\dfrac{\ln(1+a/x)}{1/x}$, form $0/0$. Differentiate: $\dfrac{-a/x^2\cdot 1/(1+a/x)}{-1/x^2}=\dfrac{a}{1+a/x}\to a$. So the limit is $e^a$.
>
> **Checks.** (a) at $x=10^{-3}$: $0.500167$; (b) at $x=50$: $2.4\times10^{-17}$; (c) at $x=10^{-6}$: $-1.4\times10^{-5}$; (d) $a=1$, $x=10^6$: $2.7182805$. ✓

> [!example] Worked Example 4.11 — The critically damped limit
> Two first-order stages in cascade with time constants $1/a$ and $1/b$ give an impulse response proportional to
> $$
> h(t)=\frac{e^{-at}-e^{-bt}}{b-a}\qquad(a\neq b).
> $$
> What happens as $a\to b$ (equal poles)?
>
> For fixed $t$, treat this as a function of $a$. Numerator and denominator both tend to $0$. Differentiate with respect to $a$: numerator derivative $-te^{-at}$, denominator derivative $-1$. So
> $$
> \lim_{a\to b}h(t)=\lim_{a\to b}\frac{-te^{-at}}{-1}=t\,e^{-bt}.
> $$
> This explains the appearance of the factor $t$ in the solution of an ODE with a repeated characteristic root ([Ch. 11](../part-4-differential-equations/11-ordinary-differential-equations.md)): it is the limit of two nearly equal exponentials.
>
> **Check.** $b=1$, $a=1.001$, $t=2$: $(e^{-2.002}-e^{-2})/(-0.001)=0.27040$; $te^{-t}=0.27067$ (the difference is $O(a-b)$). ✓

> [!warning] Common Pitfall — applying L'Hôpital where it does not apply
> L'Hôpital's rule is only for $0/0$ or $\infty/\infty$. $\lim_{x\to0}\frac{\cos x}{x+1}=1$ by substitution, but "differentiating" gives $\frac{-\sin x}{1}\to0$, which is wrong. Also, the rule differentiates numerator and denominator **separately**; it is not the quotient rule. Finally, if $\lim f'/g'$ fails to exist the rule says nothing: $\lim_{x\to\infty}(x+\sin x)/x=1$, although $(1+\cos x)/1$ has no limit.

> [!tip] Engineering Note — derivatives from measured data
> Instruments give samples, not formulas. A measured current $i[n]=i(nT_s)$ is differentiated numerically by $\big(i[n+1]-i[n-1]\big)/(2T_s)$, the central difference, which is the definition of the derivative with $h=T_s$ (made symmetric for accuracy; see §4.7). Differentiation amplifies high-frequency noise, because $j\omega$ grows with $\omega$: a $1\,$mV noise component at $1\,$MHz contributes $2\pi\times10^6\times10^{-3}\approx6.3\times10^3\,$V/s to $dv/dt$. Practical differentiators therefore band-limit the signal first.

## 4.3 Chain rule

### Motivation

Most functions in engineering are compositions. A decaying sinusoid $e^{-t/\tau}\cos(\omega t+\phi)$ has the time variable buried inside an exponential and a cosine. The gain of a filter in decibels is a logarithm of a square root of a rational function of $\omega$. A diode current is an exponential of a voltage that is itself a function of time. The chain rule differentiates compositions, and it is the most frequently used rule in calculus.

### Statement

> [!theorem] Chain rule
> If $g$ is differentiable at $x$ and $f$ is differentiable at $u=g(x)$, then $f\circ g$ is differentiable at $x$ and
> $$
> \frac{d}{dx}f\big(g(x)\big)=f'\big(g(x)\big)\,g'(x).
> $$
> In Leibniz notation, with $y=f(u)$ and $u=g(x)$:
> $$
> \frac{dy}{dx}=\frac{dy}{du}\,\frac{du}{dx}.
> $$

The Leibniz form looks like cancelling $du$, and that is a reliable mnemonic. Rates multiply along a chain: if $y$ changes 3 times as fast as $u$, and $u$ changes 2 times as fast as $x$, then $y$ changes 6 times as fast as $x$. A cascade of amplifiers with small-signal gains $A_1$ and $A_2$ has overall gain $A_1A_2$ for exactly this reason.

*Proof sketch.* The naive argument writes

$$
\frac{f(g(x+h))-f(g(x))}{h}=\frac{f(g(x+h))-f(g(x))}{g(x+h)-g(x)}\cdot\frac{g(x+h)-g(x)}{h},
$$

and lets $h\to0$. This fails if $g(x+h)=g(x)$ for small $h\neq 0$ (division by zero). The fix is to define

$$
\phi(k)=\begin{cases}\dfrac{f(u+k)-f(u)}{k}, & k\neq0\\[2mm] f'(u), & k=0\end{cases}
$$

which is continuous at $k=0$ because $f$ is differentiable at $u$. Then $f(u+k)-f(u)=\phi(k)\,k$ holds for *all* $k$, including $k=0$. Put $k=g(x+h)-g(x)$:

$$
\frac{f(g(x+h))-f(g(x))}{h}=\phi(k)\,\frac{g(x+h)-g(x)}{h}\to f'(u)\,g'(x),
$$

since $k\to0$ as $h\to0$ (because $g$ is continuous). $\blacksquare$

### Standard consequences

Applying the chain rule to each entry in the table of §4.2 gives the "generalized" rules, with $u=u(x)$:

$$
\begin{aligned}
\frac{d}{dx}u^n&=nu^{n-1}u', & \frac{d}{dx}e^{u}&=e^{u}u', & \frac{d}{dx}\ln|u|&=\frac{u'}{u},\\
\frac{d}{dx}\sin u&=\cos u\;u', & \frac{d}{dx}\cos u&=-\sin u\;u', & \frac{d}{dx}\arctan u&=\frac{u'}{1+u^2}.
\end{aligned}
$$

For a sinusoid with angular frequency $\omega$ and phase $\phi$:

$$
\frac{d}{dt}A\sin(\omega t+\phi)=A\omega\cos(\omega t+\phi)=A\omega\sin\!\left(\omega t+\phi+\frac{\pi}{2}\right).
$$

Differentiating a sinusoid multiplies its amplitude by $\omega$ and advances its phase by $90^\circ$. That is the time-domain meaning of the factor $j\omega=\omega\angle 90^\circ$. For an inductor, the voltage leads the current by $90^\circ$; for a capacitor, $i=C\,dv/dt$ means the current leads the voltage by $90^\circ$.

> [!example] Worked Example 4.12 — A damped sinusoid
> The voltage across an underdamped RLC circuit is $v(t)=V_0e^{-t/\tau}\cos(\omega_d t)$ with $V_0=10\,$V, $\tau=1\,$ms, $\omega_d=2\pi\times1000\,$rad/s. Find $dv/dt$ and evaluate it at $t=0$.
>
> **Product rule, then chain rule on each factor.** $\frac{d}{dt}e^{-t/\tau}=-\frac1\tau e^{-t/\tau}$ and $\frac{d}{dt}\cos\omega_dt=-\omega_d\sin\omega_dt$. So
> $$
> \frac{dv}{dt}=-V_0e^{-t/\tau}\left[\frac{1}{\tau}\cos\omega_dt+\omega_d\sin\omega_dt\right].
> $$
>
> **At $t=0$.** $\dfrac{dv}{dt}(0)=-V_0/\tau=-10/10^{-3}=-10^4\,$V/s $=-10\,$mV/µs.
>
> **Compact form.** With $\tan\psi=\omega_d\tau$, the bracket is $\sqrt{1/\tau^2+\omega_d^2}\,\cos(\omega_dt-\psi)$. The derivative is another damped sinusoid with the same decay, scaled by $\sqrt{\tau^{-2}+\omega_d^2}=6362\,\text{s}^{-1}$ and phase-shifted. In complex form this is immediate: $v=\operatorname{Re}\{V_0e^{st}\}$ with $s=-1/\tau+j\omega_d$, so $\dot v=\operatorname{Re}\{sV_0e^{st}\}$ and $|s|=6362\,\text{s}^{-1}$.

### Logarithmic differentiation

For products, quotients and powers of many factors, take logarithms first. If $y=f_1^{a_1}f_2^{a_2}\cdots f_n^{a_n}$ then $\ln|y|=\sum a_k\ln|f_k|$, and differentiating with the chain rule:

$$
\frac{y'}{y}=a_1\frac{f_1'}{f_1}+a_2\frac{f_2'}{f_2}+\cdots+a_n\frac{f_n'}{f_n}.
$$

The quantity $y'/y=(\ln y)'$ is the **logarithmic derivative**: the fractional (relative) rate of change. It turns products into sums and is the engine of the relative sensitivity in §4.11.

**The general power rule.** For $y=x^r$ with real $r$ and $x>0$: $\ln y=r\ln x$, so $y'/y=r/x$ and $y'=rx^{r-1}$. This completes the proof of the power rule for all real exponents.

**Variable exponents.** $y=x^x$ ($x>0$): $\ln y=x\ln x$, so $y'/y=\ln x+1$ and $y'=x^x(1+\ln x)$. Similarly $\frac{d}{dx}a^x=a^x\ln a$, from $a^x=e^{x\ln a}$.

> [!example] Worked Example 4.13 — Slope of a Bode magnitude plot
> For the first-order low-pass filter $|H|=\big(1+(\omega/\omega_c)^2\big)^{-1/2}$, the gain in decibels is $G=20\log_{10}|H|=-10\log_{10}\!\big(1+(\omega/\omega_c)^2\big)$. A Bode plot uses $x=\log_{10}\omega$ on the horizontal axis. Find the slope $dG/dx$ in dB per decade.
>
> **Chain rule.** $dG/dx=(dG/d\omega)(d\omega/dx)$. Since $\omega=10^x$, $d\omega/dx=10^x\ln10=\omega\ln10$. And
> $$
> \frac{dG}{d\omega}=-\frac{10}{\ln10}\cdot\frac{2\omega/\omega_c^2}{1+(\omega/\omega_c)^2}.
> $$
> Multiplying,
> $$
> \frac{dG}{dx}=-20\,\frac{(\omega/\omega_c)^2}{1+(\omega/\omega_c)^2}\quad\text{dB/decade}.
> $$
>
> **Interpretation.** For $\omega\ll\omega_c$ the slope is $\approx0$; at $\omega=\omega_c$ it is exactly $-10$ dB/decade; for $\omega\gg\omega_c$ it tends to $-20$ dB/decade. The asymptotic straight-line Bode plot replaces this smooth curve by its two limits, which is accurate except near the corner.
>
> **Check.** $\omega=10\omega_c$: $-20\times100/101=-19.80$ dB/decade. Numerically, $G(10.01\omega_c)-G(9.99\omega_c)$ divided by $\log_{10}(10.01/9.99)$ gives $-19.80$. ✓

> [!tip] Engineering Note — why decibels and log axes
> The logarithmic derivative explains the habit of plotting gain in dB against $\log\omega$. If $|H|\propto\omega^{n}$ over some band, then $\log|H|=n\log\omega+\text{const}$, a straight line of slope $n$ (i.e. $20n$ dB/decade). Each pole contributes $-20$ dB/decade and each zero $+20$ dB/decade beyond its corner frequency because logarithms turn the product of factors in $H$ into a sum of their contributions. The slope $d\ln|H|/d\ln\omega$ is exactly the relative sensitivity $S^{|H|}_{\omega}$ of §4.11.

### Derivatives of inverse functions

If $y=f(x)$ is differentiable and strictly monotonic with $f'(x)\ne0$, its inverse $x=f^{-1}(y)$ is differentiable and

$$
\frac{d}{dy}f^{-1}(y)=\frac{1}{f'(x)}\bigg|_{x=f^{-1}(y)},\qquad\text{i.e.}\qquad\frac{dx}{dy}=\frac{1}{dy/dx}.
$$

*Derivation.* Differentiate the identity $f\big(f^{-1}(y)\big)=y$ with respect to $y$ using the chain rule: $f'\big(f^{-1}(y)\big)\cdot(f^{-1})'(y)=1$. Graphically, the graph of $f^{-1}$ is the mirror image of the graph of $f$ in the line $y=x$, which turns a slope $m$ into $1/m$.

**Arctangent.** $y=\arctan x$ means $x=\tan y$ with $y\in(-\pi/2,\pi/2)$. Then $dx/dy=\sec^2y=1+\tan^2y=1+x^2$, so $dy/dx=1/(1+x^2)$.

**Arcsine.** $x=\sin y$, $y\in[-\pi/2,\pi/2]$: $dx/dy=\cos y=+\sqrt{1-x^2}$ (positive on that interval), so $dy/dx=1/\sqrt{1-x^2}$.

**Logarithm.** $x=e^y$: $dx/dy=e^y=x$, so $\frac{d}{dx}\ln x=1/x$ (a second proof).

> [!example] Worked Example 4.14 — Inverse of a device characteristic
> A nonlinear sensor has output $y=f(x)=x^3+x$ (in normalized units). It is used "in reverse": the measured $y$ is converted to $x=f^{-1}(y)$. How sensitive is the recovered $x$ to errors in $y$ when $y=2$?
>
> $f(1)=2$, so $f^{-1}(2)=1$. $f'(x)=3x^2+1$, so $f'(1)=4$. Hence $(f^{-1})'(2)=1/4$: an error of $0.01$ in $y$ produces an error of only about $0.0025$ in $x$.
>
> **Check.** Solve $x^3+x=2.01$ numerically: $x=1.0024953$; the change is $0.0024953\approx0.01/4$. ✓
>
> We never needed a formula for $f^{-1}$, which would require solving a cubic.

### Parametric curves

If a curve is given as $x=x(t)$, $y=y(t)$ (a phasor tip, a Lissajous figure on an oscilloscope, a trajectory in a phase plane), then wherever $dx/dt\neq0$,

$$
\frac{dy}{dx}=\frac{dy/dt}{dx/dt}.
$$

This is the chain rule $\frac{dy}{dt}=\frac{dy}{dx}\frac{dx}{dt}$ solved for $dy/dx$. The second derivative needs one more application:

$$
\frac{d^2y}{dx^2}=\frac{d}{dx}\left(\frac{dy}{dx}\right)=\frac{\frac{d}{dt}\left(\frac{dy}{dx}\right)}{dx/dt}.
$$

> [!example] Worked Example 4.15 — Slope of a Lissajous figure
> An oscilloscope in X–Y mode displays $x=\cos t$, $y=\sin 2t$. Find the slope of the trace at $t=\pi/6$ and the points where the tangent is horizontal or vertical.
>
> $dx/dt=-\sin t$, $dy/dt=2\cos2t$. At $t=\pi/6$: $dx/dt=-\tfrac12$, $dy/dt=2\cdot\tfrac12=1$. Slope $=1/(-\tfrac12)=-2$, at the point $(\cos\frac\pi6,\sin\frac\pi3)=(0.866,0.866)$.
>
> **Horizontal tangents** need $dy/dt=0$ with $dx/dt\neq0$: $\cos2t=0\Rightarrow t=\pi/4,3\pi/4,\ldots$, giving $(\pm0.707,\pm1)$, the four peaks of the figure-eight.
>
> **Vertical tangents** need $dx/dt=0$ with $dy/dt\neq0$: $t=0,\pi$, where $dy/dt=2$. The tangent is vertical at $(\pm1,0)$.

```desmos-graph
left=-1.5; right=1.5; bottom=-1.3; top=1.3
---
(\cos(t),\sin(2t))|blue
y=-2(x-0.866)+0.866|red|0.4<x<1.3
(0.866,0.866)|label:slope -2
```
*The Lissajous figure $x=\cos t$, $y=\sin 2t$ (a 1:2 frequency ratio) and its tangent at $t=\pi/6$.*

> [!example] Worked Example 4.16 — Small-signal conductance of a diode
> The Shockley diode equation is $I=I_s\big(e^{V/(nV_T)}-1\big)$. Find the incremental conductance $g_d=dI/dV$ and evaluate it at $I=1\,$mA with $nV_T=26\,$mV.
>
> **Chain rule** on $e^{V/(nV_T)}$ (inner derivative $1/(nV_T)$):
> $$
> g_d=\frac{dI}{dV}=\frac{I_s}{nV_T}e^{V/(nV_T)}=\frac{I+I_s}{nV_T}\approx\frac{I}{nV_T}.
> $$
> The last step uses $I_s\ll I$ (typically $I_s\sim10^{-14}$ A). Note that $g_d$ depends only on the bias current, not on $I_s$: a universal result.
>
> **Numbers.** $g_d=10^{-3}/0.026=38.5\,$mS, so the small-signal resistance is $r_d=1/g_d=26\,\Omega$.
>
> **If the voltage varies in time,** the chain rule again gives $\dfrac{dI}{dt}=g_d\,\dfrac{dV}{dt}$: the diode responds to small, slow changes like a resistor $r_d$. We turn this into a full linear model in §4.8.

> [!example] Worked Example 4.17 — Rate of change of an impedance magnitude
> A series RLC circuit has $|Z(\omega)|=\sqrt{R^2+X^2}$ with reactance $X(\omega)=\omega L-\dfrac{1}{\omega C}$. Find $d|Z|/d\omega$ and show that $|Z|$ has a stationary point at resonance.
>
> **Outer function** $\sqrt{\cdot}$, inner function $R^2+X^2$, and $X$ itself depends on $\omega$:
> $$
> \frac{d|Z|}{d\omega}=\frac{1}{2\sqrt{R^2+X^2}}\cdot 2X\frac{dX}{d\omega}=\frac{X}{|Z|}\left(L+\frac{1}{\omega^2C}\right).
> $$
> The bracket is always positive, so the sign of $d|Z|/d\omega$ is the sign of $X$. At resonance $\omega_0=1/\sqrt{LC}$, $X=0$, so $d|Z|/d\omega=0$: $|Z|$ has its minimum value $R$ there.
>
> **Numbers.** $R=10\,\Omega$, $L=1\,$mH, $C=1\,\mu$F: $\omega_0=31\,623\,$rad/s. At $\omega=1.1\omega_0$: $X=34.785-28.748=6.037\,\Omega$, $|Z|=11.68\,\Omega$, and $d|Z|/d\omega=(6.037/11.68)(10^{-3}+0.826\times10^{-3})=9.44\times10^{-4}\,\Omega\,$s/rad.

> [!warning] Common Pitfall — forgetting the inner derivative, and degrees
> $\frac{d}{dx}\sin(x^2)=2x\cos(x^2)$, not $\cos(x^2)$ and not $\cos(2x)$. Every layer of a composition contributes a factor. A subtler instance: $\frac{d}{dx}\sin(x^\circ)$, with $x$ in degrees, is $\frac{\pi}{180}\cos(x^\circ)$, because $x^\circ=\pi x/180$ radians. Calculators set to degrees silently break every derivative formula in this chapter.

## 4.4 Implicit differentiation

### Motivation

Many relations in engineering cannot be solved explicitly for the quantity of interest. The current in a resistor–diode circuit satisfies

$$
V_s=IR+nV_T\ln\!\left(1+\frac{I}{I_s}\right),
$$

which has no closed-form solution for $I$ in elementary functions. Yet we still want $dI/dV_s$: how much the current moves when the supply moves. Implicit differentiation finds derivatives directly from an equation $F(x,y)=0$ without solving it.

### Method

Treat $y$ as an (unknown) differentiable function of $x$, differentiate **both sides** of the equation with respect to $x$, applying the chain rule to every term containing $y$ (so $\frac{d}{dx}y^2=2y\,y'$, $\frac{d}{dx}\sin y=\cos y\,y'$, $\frac{d}{dx}(xy)=y+xy'$), then solve the resulting linear equation for $y'$.

The result is always linear in $y'$, because the chain rule produces exactly one factor $y'$ per occurrence of $y$.

> [!example] Worked Example 4.18 — A circle
> Find $dy/dx$ and $d^2y/dx^2$ on $x^2+y^2=r^2$.
>
> **First derivative.** $2x+2y\,y'=0\Rightarrow y'=-x/y$ (valid where $y\neq0$). At $(r/\sqrt2,r/\sqrt2)$ the slope is $-1$. The tangent is perpendicular to the radius, whose slope is $y/x$; the product of slopes is $-1$. ✓
>
> **Second derivative.** Differentiate $y'=-x/y$ by the quotient rule, then substitute for $y'$:
> $$
> y''=-\frac{y-xy'}{y^2}=-\frac{y+x^2/y}{y^2}=-\frac{x^2+y^2}{y^3}=-\frac{r^2}{y^3}.
> $$
> On the upper half ($y>0$) the circle is concave down ($y''<0$), as it should be.
>
> **Check against the explicit form** $y=\sqrt{r^2-x^2}$: $y'=-x/\sqrt{r^2-x^2}=-x/y$. ✓

> [!example] Worked Example 4.19 — Tangent to the folium of Descartes
> Find the tangent line to $x^3+y^3=6xy$ at $(3,3)$.
>
> First check the point: $27+27=54=6\cdot9$. ✓ Differentiate:
> $$
> 3x^2+3y^2y'=6y+6xy'\quad\Longrightarrow\quad y'(3y^2-6x)=6y-3x^2\quad\Longrightarrow\quad y'=\frac{2y-x^2}{y^2-2x}.
> $$
> At $(3,3)$: $y'=(6-9)/(9-6)=-1$. Tangent: $y-3=-(x-3)$, i.e. $y=6-x$.
>
> **Where is $y'$ undefined?** Where $y^2=2x$. Combined with the curve equation this gives the origin and the point $(2^{5/3},2^{4/3})\approx(3.17,2.52)$, where the tangent is vertical.

```desmos-graph
left=-4; right=6; bottom=-4; top=6
---
x^3+y^3=6xy|blue
y=6-x|red|dashed
(3,3)|label:(3,3), slope -1
```
*The folium of Descartes $x^3+y^3=6xy$ is not the graph of a single function $y(x)$, yet implicit differentiation gives its slope at every regular point; the tangent at $(3,3)$ is $y=6-x$.*

### The general formula and when it is valid

For a relation $F(x,y)=0$, differentiating with respect to $x$ (using the multivariable chain rule of §4.5) gives $F_x+F_y\,y'=0$, so

$$
\frac{dy}{dx}=-\frac{F_x}{F_y}=-\frac{\partial F/\partial x}{\partial F/\partial y}\qquad(F_y\neq0).
$$

> [!theorem] Implicit function theorem (one-variable form, stated)
> If $F$ has continuous partial derivatives near $(a,b)$, $F(a,b)=0$ and $F_y(a,b)\neq0$, then near $x=a$ the equation $F(x,y)=0$ defines a unique differentiable function $y=y(x)$ with $y(a)=b$, and $y'=-F_x/F_y$.

The condition $F_y\neq0$ matters. At points where $F_y=0$ (the leftmost and rightmost points of a circle, the self-intersection of the folium at the origin), the curve may have a vertical tangent, a branch point, or may not be a function of $x$ at all. [Ch. 6](06-multivariable-calculus.md) proves the theorem and extends it to systems of equations via the Jacobian.

> [!example] Worked Example 4.20 — Resistor–diode circuit
> A source $V_s$ drives a resistor $R=1\,$kΩ in series with a diode ($nV_T=26\,$mV, $I_s=10^{-14}\,$A). KVL gives
> $$
> V_s=IR+nV_T\ln\!\left(1+\frac{I}{I_s}\right).
> $$
> Find $dI/dV_s$ at the operating point $I=1\,$mA.
>
> **Differentiate both sides with respect to $V_s$**, treating $I=I(V_s)$:
> $$
> 1=R\frac{dI}{dV_s}+nV_T\cdot\frac{1/I_s}{1+I/I_s}\frac{dI}{dV_s}=\left(R+\frac{nV_T}{I_s+I}\right)\frac{dI}{dV_s}.
> $$
> Hence
> $$
> \frac{dI}{dV_s}=\frac{1}{R+r_d},\qquad r_d=\frac{nV_T}{I+I_s}\approx\frac{nV_T}{I}.
> $$
>
> **Numbers.** $r_d=0.026/10^{-3}=26\,\Omega$, so $dI/dV_s=1/1026\,\Omega=0.975\,$mS. A $10\,$mV ripple on the supply produces about $9.75\,\mu$A of current ripple.
>
> **Interpretation.** The diode behaves, for small changes, like a $26\,\Omega$ resistor in series with $R$: exactly the small-signal model of Example 4.16, obtained here without ever solving the nonlinear equation.
>
> **Numerical check.** Solving the KVL equation numerically at $V_s\pm1\,$mV around the operating point and forming the central difference gives $0.9747\,$mS. ✓

> [!example] Worked Example 4.21 — Horizontal tangents of an ellipse
> The rotated ellipse $x^2+xy+y^2=7$ passes through $(1,2)$. Find the slope there and the points with horizontal tangents.
>
> Differentiate: $2x+(y+xy')+2yy'=0\Rightarrow y'=-\dfrac{2x+y}{x+2y}$. At $(1,2)$: $y'=-4/5$.
>
> **Horizontal tangents:** $2x+y=0$, i.e. $y=-2x$. Substitute: $x^2-2x^2+4x^2=3x^2=7$, so $x=\pm\sqrt{7/3}=\pm1.528$, $y=\mp3.055$. The ellipse's highest point is $(-1.528,3.055)$ and its lowest $(1.528,-3.055)$.
>
> **Check.** $1.528^2+1.528\cdot(-3.055)+3.055^2=2.333-4.667+9.333=7.000$. ✓

> [!tip] Engineering Note — implicit relations are everywhere in circuits
> Any circuit containing a nonlinear element, solved by nodal analysis, produces an implicit equation (or system of equations) for the node voltages. Circuit simulators such as SPICE do not solve these explicitly. They linearize them (§4.8) and iterate with Newton's method, and the derivatives they need are exactly the implicit derivatives of this section, assembled into a Jacobian matrix ([Ch. 6](06-multivariable-calculus.md), [Ch. 22](../part-7-numerical-computational/22-numerical-methods.md)).

> [!warning] Common Pitfall — differentiating $y$ as if it were constant
> When differentiating implicitly with respect to $x$, $\frac{d}{dx}(y^3)$ is $3y^2y'$, not $3y^2$ and not $0$. Likewise, $\frac{d}{dx}(xy)$ needs the product rule: $y+xy'$. Forgetting the $y'$ factor is the most common error in this topic and always produces an equation that cannot be solved for $y'$ correctly.

## 4.5 Partial derivatives

### Motivation

Engineering quantities rarely depend on one variable. The power dissipated in a resistor, $P=V^2/R$, depends on both voltage and resistance. The cutoff frequency of an RC filter, $f_c=1/(2\pi RC)$, depends on two components. A potential $\varphi(x,y,z)$ depends on three coordinates, and a travelling wave $u(x,t)$ on position and time. To ask "how does $P$ change if only $R$ drifts?" we hold the other variables fixed and differentiate with respect to one. That is a partial derivative. This section is an introduction, enough for the rest of this chapter (implicit differentiation, linearization in several variables, sensitivity analysis). [Ch. 6](06-multivariable-calculus.md) develops multivariable calculus fully: total differentials, gradients, directional derivatives, Jacobians and differentiability.

### Definition

> [!definition] Partial derivative
> Let $f(x,y)$ be defined near $(a,b)$. The **partial derivatives** of $f$ at $(a,b)$ are
> $$
> \frac{\partial f}{\partial x}(a,b)=\lim_{h\to0}\frac{f(a+h,b)-f(a,b)}{h},\qquad
> \frac{\partial f}{\partial y}(a,b)=\lim_{h\to0}\frac{f(a,b+h)-f(a,b)}{h},
> $$
> when the limits exist. Other notations: $f_x$, $\partial_x f$, $D_1f$; and $\left(\dfrac{\partial f}{\partial x}\right)_{y}$ when it is important to say which variable is held fixed. The definition extends to any number of variables: vary one, freeze all the others.

The symbol $\partial$ ("partial", "del", or "dee") distinguishes a partial derivative from an ordinary one. Computing $\partial f/\partial x$ needs no new rules: treat every other variable as a constant and use the one-variable rules of §4.2–4.3.

**Geometric meaning.** The graph of $f(x,y)$ is a surface $z=f(x,y)$. Freezing $y=b$ slices the surface with the vertical plane $y=b$, leaving a curve $z=f(x,b)$. Then $f_x(a,b)$ is the slope of that curve at $x=a$: the slope of the surface in the $x$-direction. Similarly, $f_y$ is the slope in the $y$-direction. A topographic analogy: $f_x$ is how steeply the ground rises as you walk due east, $f_y$ as you walk due north.

> [!example] Worked Example 4.22 — Computing partial derivatives
> Find $f_x$ and $f_y$ for (a) $f=x^2y^3+\sin(xy)$, (b) $f=e^{-x}\cos(3y)$, (c) $f=\sqrt{x^2+y^2}$.
>
> (a) Holding $y$ fixed, $x^2y^3$ is ($y^3$ times) $x^2$, and $\sin(xy)$ needs the chain rule with inner derivative $\partial(xy)/\partial x=y$:
> $$
> f_x=2xy^3+y\cos(xy),\qquad f_y=3x^2y^2+x\cos(xy).
> $$
>
> (b) $f_x=-e^{-x}\cos3y$, $f_y=-3e^{-x}\sin3y$.
>
> (c) With $r=\sqrt{x^2+y^2}$: $f_x=\dfrac{x}{\sqrt{x^2+y^2}}=\dfrac xr$, $f_y=\dfrac yr$ (for $r\neq0$). At $(3,4)$: $f_x=0.6$, $f_y=0.8$. These are the direction cosines of the point: the distance from the origin grows fastest in the radial direction.

> [!example] Worked Example 4.23 — Partial derivatives of circuit formulas
> (a) For $P=V^2/R$: $\dfrac{\partial P}{\partial V}=\dfrac{2V}{R}$ (units W/V $=$ A) and $\dfrac{\partial P}{\partial R}=-\dfrac{V^2}{R^2}$ (W/Ω). At $V=12\,$V, $R=6\,\Omega$: $\partial P/\partial V=4\,$A and $\partial P/\partial R=-4\,$W/Ω.
>
> (b) Two resistors in parallel, $R_p=\dfrac{R_1R_2}{R_1+R_2}$. By the quotient rule with $R_2$ constant:
> $$
> \frac{\partial R_p}{\partial R_1}=\frac{R_2(R_1+R_2)-R_1R_2}{(R_1+R_2)^2}=\frac{R_2^2}{(R_1+R_2)^2},\qquad
> \frac{\partial R_p}{\partial R_2}=\frac{R_1^2}{(R_1+R_2)^2}.
> $$
> With $R_1=1\,$kΩ, $R_2=4\,$kΩ ($R_p=800\,\Omega$): $\partial R_p/\partial R_1=16/25=0.64$ and $\partial R_p/\partial R_2=1/25=0.04$. The parallel combination is dominated by the smaller resistor; a $10\,\Omega$ change in $R_1$ moves $R_p$ by about $6.4\,\Omega$, the same change in $R_2$ by only $0.4\,\Omega$.
>
> **Check.** $R_1=1010$: $R_p=1010\cdot4000/5010=806.39\,\Omega$, an increase of $6.39\,\Omega$. ✓

### The first-order change: a preview of the total differential

If both $V$ and $R$ change a little, the change in $P$ is approximately the sum of the separate effects:

$$
\Delta P\approx\frac{\partial P}{\partial V}\Delta V+\frac{\partial P}{\partial R}\Delta R .
$$

*Why.* Change $V$ first and then $R$:

$$
\Delta P=\big[P(V+\Delta V,R+\Delta R)-P(V,R+\Delta R)\big]+\big[P(V,R+\Delta R)-P(V,R)\big].
$$

By the mean value theorem, each bracket is a one-variable derivative times the step, evaluated at an intermediate point. If the partials are continuous, the intermediate points can be replaced by $(V,R)$ at the cost of an error that is small compared with $|\Delta V|+|\Delta R|$.

With $V=12\,$V, $R=6\,\Omega$, $\Delta V=+0.2\,$V and $\Delta R=-0.1\,\Omega$: $\Delta P\approx4(0.2)+(-4)(-0.1)=1.2\,$W. The exact change is $12.2^2/5.9-24=1.227\,$W. The small discrepancy comes from second-order terms. The expression $df=f_x\,dx+f_y\,dy$ is the **total differential**; it underlies linearization (§4.8) and sensitivity analysis (§4.11), and is studied properly in [Ch. 6](06-multivariable-calculus.md).

### The chain rule for several variables

If $z=f(x,y)$ and both $x=x(t)$ and $y=y(t)$ depend on $t$, dividing the first-order change by $\Delta t$ and letting $\Delta t\to0$ gives

$$
\frac{dz}{dt}=\frac{\partial f}{\partial x}\frac{dx}{dt}+\frac{\partial f}{\partial y}\frac{dy}{dt}.
$$

Each path from $t$ to $z$ contributes a product of derivatives, and the contributions add. Applied to $F(x,y(x))=0$ with $t=x$, this gives $F_x+F_y\,y'=0$, the implicit-differentiation formula of §4.4.

> [!example] Worked Example 4.24 — A heating resistor on a drifting supply
> A resistor is fed from a supply whose voltage is rising at $dV/dt=0.1\,$V/s while the resistor heats and its resistance rises at $dR/dt=0.02\,\Omega$/s. At the instant when $V=12\,$V and $R=6\,\Omega$, how fast is the dissipated power changing?
>
> $$
> \frac{dP}{dt}=\frac{\partial P}{\partial V}\frac{dV}{dt}+\frac{\partial P}{\partial R}\frac{dR}{dt}=\frac{2V}{R}(0.1)-\frac{V^2}{R^2}(0.02)=4(0.1)-4(0.02)=0.32\ \text{W/s}.
> $$
>
> The voltage rise increases the power by $0.40\,$W/s; the resistance rise reduces it by $0.08\,$W/s.
>
> **Units.** A·(V/s) $=$ W/s and (W/Ω)(Ω/s) $=$ W/s. ✓

### Higher partial derivatives and Clairaut's theorem

Partial derivatives can be differentiated again. With two variables there are four second partials:

$$
f_{xx}=\frac{\partial^2 f}{\partial x^2},\qquad f_{yy}=\frac{\partial^2 f}{\partial y^2},\qquad f_{xy}=\frac{\partial}{\partial y}\left(\frac{\partial f}{\partial x}\right)=\frac{\partial^2 f}{\partial y\,\partial x},\qquad f_{yx}=\frac{\partial^2 f}{\partial x\,\partial y}.
$$

> [!theorem] Clairaut's (Schwarz's) theorem
> If $f_{xy}$ and $f_{yx}$ are both continuous near $(a,b)$, then $f_{xy}(a,b)=f_{yx}(a,b)$. The order of differentiation does not matter.

For the function of Example 4.22(a): $f_{xy}=\partial_y\big(2xy^3+y\cos(xy)\big)=6xy^2+\cos(xy)-xy\sin(xy)$, and $f_{yx}=\partial_x\big(3x^2y^2+x\cos(xy)\big)=6xy^2+\cos(xy)-xy\sin(xy)$. They agree.

Clairaut's theorem is behind several facts in field theory. For example, the curl of a gradient vanishes, $\nabla\times\nabla\varphi=\mathbf{0}$, because each component is a difference of equal mixed partials ([Ch. 7](07-vector-calculus.md)).

> [!example] Worked Example 4.25 — Checking that a function solves a PDE
> (a) Show that $\varphi=\ln(x^2+y^2)$ satisfies Laplace's equation $\varphi_{xx}+\varphi_{yy}=0$ for $(x,y)\neq(0,0)$.
>
> $\varphi_x=\dfrac{2x}{x^2+y^2}$. Quotient rule: $\varphi_{xx}=\dfrac{2(x^2+y^2)-2x\cdot2x}{(x^2+y^2)^2}=\dfrac{2(y^2-x^2)}{(x^2+y^2)^2}$. By the symmetry $x\leftrightarrow y$, $\varphi_{yy}=\dfrac{2(x^2-y^2)}{(x^2+y^2)^2}$. The sum is $0$. ✓ (Up to a constant factor this is the potential of an infinite line charge.)
>
> (b) Show that $u(x,t)=\sin(x-ct)$ satisfies the wave equation $u_{tt}=c^2u_{xx}$.
>
> $u_x=\cos(x-ct)$, $u_{xx}=-\sin(x-ct)$; $u_t=-c\cos(x-ct)$, $u_{tt}=-c^2\sin(x-ct)$. So $u_{tt}=c^2u_{xx}$. ✓ The same holds for *any* twice-differentiable $g(x-ct)$, by the chain rule: a waveform travelling at speed $c$.
>
> (c) Show that $u=e^{-\alpha k^2t}\sin(kx)$ satisfies the heat (diffusion) equation $u_t=\alpha u_{xx}$.
>
> $u_t=-\alpha k^2e^{-\alpha k^2t}\sin kx$ and $u_{xx}=-k^2e^{-\alpha k^2t}\sin kx$, so $u_t=\alpha u_{xx}$. ✓ High spatial frequencies ($k$ large) decay fastest, as $e^{-\alpha k^2 t}$.
>
> These checks are the starting point for [Ch. 13](../part-4-differential-equations/13-partial-differential-equations.md).

> [!tip] Engineering Note — fields from potentials
> The electric field is the negative gradient of the potential: $E_x=-\partial\varphi/\partial x$, $E_y=-\partial\varphi/\partial y$, $E_z=-\partial\varphi/\partial z$. For a point charge, $\varphi=\dfrac{q}{4\pi\varepsilon_0}(x^2+y^2+z^2)^{-1/2}$, and the chain rule gives $E_x=\dfrac{q}{4\pi\varepsilon_0}\dfrac{x}{R^3}$ with $R=\sqrt{x^2+y^2+z^2}$, and similarly for $E_y$, $E_z$. Together these give the familiar inverse-square law $|\mathbf{E}|=q/(4\pi\varepsilon_0R^2)$, pointing radially. Partial derivatives turn one scalar function into the three components of a vector field.

> [!warning] Common Pitfall — partial derivatives existing is not enough
> In one variable, differentiable implies continuous. In two variables, the mere existence of $f_x$ and $f_y$ does **not** imply continuity. The function $f(x,y)=\dfrac{xy}{x^2+y^2}$ (with $f(0,0)=0$) has $f_x(0,0)=f_y(0,0)=0$, because $f$ is identically zero on both axes. Yet along the line $y=x$, $f=\tfrac12$ everywhere except the origin, so $f$ is not continuous at $(0,0)$. The proper notion of differentiability in several variables, and the condition (continuous partials) that guarantees it, is in [Ch. 6](06-multivariable-calculus.md). Also be explicit about which variables are held fixed: in thermodynamics and device physics, $(\partial f/\partial x)_y$ and $(\partial f/\partial x)_z$ are different quantities.

## 4.6 Higher-order derivatives

### Definitions and meaning

The derivative $f'$ is itself a function and can be differentiated again. The **second derivative** is $f''=(f')'$, written $\dfrac{d^2f}{dx^2}$ or $\ddot f$ for time derivatives; the **$n$-th derivative** is $f^{(n)}=\dfrac{d^nf}{dx^n}$.

Physically, if $x(t)$ is position, $\dot x$ is velocity, $\ddot x$ acceleration and $\dddot x$ jerk. In circuits, if $q(t)$ is the charge on a capacitor, $\dot q=i$ is current and $\ddot q=di/dt$ determines the inductor voltage $L\ddot q$. The series RLC equation

$$
L\ddot q+R\dot q+\frac{1}{C}q=v_s(t)
$$

is a statement about the zeroth, first and second derivatives of one function. Every ODE of [Ch. 11](../part-4-differential-equations/11-ordinary-differential-equations.md) and Ch. 12 is such a statement.

Geometrically, $f''$ measures how the slope changes: **concavity**.

> [!definition] Concavity and inflection points
> On an interval where $f''>0$, the slope $f'$ is increasing and the graph is **concave up** (it lies above its tangent lines). Where $f''<0$, the graph is **concave down**. A point where $f$ is continuous and the concavity changes is an **inflection point**; there $f''=0$ or $f''$ fails to exist.

> [!example] Worked Example 4.26 — Shape of $te^{-t}$
> The pulse $f(t)=te^{-t}$ (the critically damped response of Example 4.11 with $b=1$) appears constantly in circuit transients. Find $f'$, $f''$, the peak and the inflection point.
>
> $$
> f'=(1-t)e^{-t},\qquad f''=-e^{-t}-(1-t)e^{-t}=(t-2)e^{-t}.
> $$
> - $f'=0$ at $t=1$; $f''(1)=-e^{-1}<0$, so this is a maximum, $f(1)=1/e=0.368$.
> - $f''=0$ at $t=2$ and changes sign there (from $-$ to $+$), so $t=2$ is an inflection point, with $f(2)=2e^{-2}=0.271$. After $t=2$ the pulse is concave up and decays towards zero.
> - At the inflection point $f'(2)=-e^{-2}=-0.135$: the pulse is falling fastest there.

```desmos-graph
left=-0.5; right=7; bottom=-1.2; top=1.2
---
y=xe^{-x}|blue|x>0
y=(1-x)e^{-x}|red|dashed|x>0
y=(x-2)e^{-x}|green|dotted|x>0
(1,0.368)|label:peak
(2,0.271)|label:inflection
```
*The pulse $te^{-t}$ (blue), its first derivative (red, dashed; zero at the peak) and its second derivative (green, dotted; zero at the inflection point, where the decay is steepest).*

### Patterns for $n$-th derivatives

| $f(x)$ | $f^{(n)}(x)$ |
|---|---|
| $x^m$ ($m$ a positive integer) | $\dfrac{m!}{(m-n)!}x^{m-n}$ for $n\le m$; $0$ for $n>m$ |
| $e^{ax}$ | $a^ne^{ax}$ |
| $\sin(ax)$ | $a^n\sin\!\left(ax+\dfrac{n\pi}{2}\right)$ |
| $\cos(ax)$ | $a^n\cos\!\left(ax+\dfrac{n\pi}{2}\right)$ |
| $\dfrac{1}{1-x}$ | $\dfrac{n!}{(1-x)^{n+1}}$ |
| $\ln x$ | $(-1)^{n-1}\dfrac{(n-1)!}{x^n}$ ($n\ge1$) |

Each row is proved by induction: compute the first few derivatives, guess the pattern, and check that differentiating the $n$-th gives the $(n+1)$-th. For example, $\frac{d}{dx}\,n!(1-x)^{-(n+1)}=n!(n+1)(1-x)^{-(n+2)}=(n+1)!(1-x)^{-(n+2)}$. The sine row says that each differentiation multiplies by $a$ and shifts the phase by $90^\circ$, the same fact as $(j\omega)^n$ in phasor form.

### Leibniz's rule for products

> [!theorem] General Leibniz rule
> If $f$ and $g$ are $n$ times differentiable,
> $$
> (fg)^{(n)}=\sum_{k=0}^{n}\binom{n}{k}f^{(k)}g^{(n-k)} .
> $$

*Proof sketch.* For $n=1$ this is the product rule. Assume it for $n$ and differentiate each term with the product rule; the coefficients combine by Pascal's identity $\binom nk+\binom n{k-1}=\binom{n+1}{k}$, exactly as in the proof of the binomial theorem. The rule is the binomial theorem with powers replaced by derivatives.

> [!example] Worked Example 4.27 — $n$-th derivative of $x^2e^{2x}$
> With $f=x^2$ (only $f,f',f''$ nonzero: $x^2,2x,2$) and $g=e^{2x}$ ($g^{(m)}=2^me^{2x}$), only three terms of Leibniz's rule survive:
> $$
> \frac{d^n}{dx^n}\big(x^2e^{2x}\big)=x^2\,2^ne^{2x}+n\,(2x)\,2^{n-1}e^{2x}+\binom n2\,2\cdot2^{n-2}e^{2x}
> =2^{n-2}e^{2x}\big(4x^2+4nx+n(n-1)\big).
> $$
>
> **Checks.** $n=1$: $2^{-1}e^{2x}(4x^2+4x)=e^{2x}(2x^2+2x)$, which matches the product rule. $n=3$: $2e^{2x}(4x^2+12x+6)=e^{2x}(8x^2+24x+12)$, which matches direct differentiation (verified symbolically for $n=1,\ldots,5$). ✓

> [!example] Worked Example 4.28 — Verifying a circuit response
> Show that $q(t)=e^{-t}\sin2t$ satisfies $\ddot q+2\dot q+5q=0$, the source-free equation of a series RLC circuit with (normalized) $L=1$, $R=2$, $1/C=5$.
>
> $$
> \begin{aligned}
> \dot q&=e^{-t}(-\sin2t+2\cos2t),\\
> \ddot q&=e^{-t}\big(\sin2t-2\cos2t\big)+e^{-t}\big(-2\cos2t-4\sin2t\big)=e^{-t}(-3\sin2t-4\cos2t).
> \end{aligned}
> $$
> Then
> $$
> \ddot q+2\dot q+5q=e^{-t}\big[(-3-2+5)\sin2t+(-4+4+0)\cos2t\big]=0. \checkmark
> $$
>
> The faster route uses complex exponentials: $q=\operatorname{Im}\{e^{st}\}$ with $s=-1+2j$, and each derivative multiplies by $s$, so the equation becomes $s^2+2s+5=0$. Indeed $(-1+2j)^2+2(-1+2j)+5=(1-4j-4)+(-2+4j)+5=0$. Turning derivatives into powers of $s$ is the central idea of [Ch. 12](../part-4-differential-equations/12-linear-differential-equations.md) and [Ch. 14](../part-5-transforms/14-laplace-transforms.md).

### Smoothness classes

A function is of class $C^k$ on an interval if its first $k$ derivatives exist and are continuous there; $C^\infty$ means all derivatives exist (polynomials, $e^x$, $\sin x$). The hierarchy is strict:

- $|x|$ is $C^0$ but not differentiable at $0$;
- $x|x|$ has derivative $2|x|$, so it is $C^1$ but not $C^2$;
- $x^2|x|$ is $C^2$ but not $C^3$.

Smoothness matters in practice: a waveform whose $k$-th derivative jumps has a Fourier spectrum that decays like $1/\omega^{k+1}$ ([Ch. 15](../part-5-transforms/15-fourier-series.md)). A square wave (a jump in the function itself) has harmonics falling as $1/n$; a triangle wave (a jump in the first derivative) as $1/n^2$. Smoother edges mean less high-frequency interference, which is why digital edge rates are deliberately controlled.

### Curvature

The second derivative alone does not measure how sharply a curve bends, because it depends on how steep the curve is. The **curvature** of $y=f(x)$ is

$$
\kappa=\frac{|f''(x)|}{\big(1+f'(x)^2\big)^{3/2}},
$$

and the radius of curvature is $1/\kappa$. For $y=x^2$ at the origin, $\kappa=2$, so the best-fitting circle there has radius $\tfrac12$. For a circle of radius $r$, Example 4.18 gave $y''=-r^2/y^3$ and $y'=-x/y$, so $\kappa=\dfrac{r^2/|y|^3}{(1+x^2/y^2)^{3/2}}=\dfrac{r^2/|y|^3}{r^3/|y|^3}=\dfrac1r$, as it must be.

> [!tip] Engineering Note — second derivatives in measurement
> The second derivative is a sensitive detector of features that are hard to see in the raw data: the "knee" of a diode or MOSFET characteristic, the onset of saturation in a magnetic core, or a peak position in a spectrum. Numerically, it is estimated by the central second difference
> $$
> f''(x)\approx\frac{f(x+h)-2f(x)+f(x-h)}{h^2},
> $$
> whose error is $O(h^2)$ (derived from Taylor series in §4.7). Because of the $h^2$ in the denominator, noise of size $\sigma$ in $f$ produces an error of order $\sigma/h^2$ in $f''$: second derivatives of measured data must be heavily smoothed.

> [!warning] Common Pitfall — $f''(c)=0$ does not guarantee an inflection point
> For $f=x^4$, $f''(0)=0$, but $f''=12x^2\ge0$ on both sides: no change of concavity, so no inflection point ($x=0$ is in fact a minimum). Always check that $f''$ changes sign. Also note that $\dfrac{d^2y}{dx^2}\neq\left(\dfrac{dy}{dx}\right)^2$ and, for parametric curves, $\dfrac{d^2y}{dx^2}\neq\dfrac{d^2y/dt^2}{d^2x/dt^2}$ (use the formula in §4.3).

## 4.7 Taylor and Maclaurin series

### Motivation

Polynomials are the only functions a processor can evaluate directly: they need nothing but additions and multiplications. Taylor's idea is to replace a complicated function, near a chosen point, by the polynomial that matches its value and as many derivatives as we like. This single idea gives:

- the small-signal models of diodes and transistors (keep the linear term),
- harmonic-distortion estimates (keep the quadratic and cubic terms),
- the approximations $\sin\theta\approx\theta$, $\cos\theta\approx1-\theta^2/2$, $(1+x)^n\approx1+nx$ used throughout engineering,
- how calculators and DSP libraries compute $e^x$, $\sin x$ and $\ln x$,
- the error analysis of every numerical method in [Ch. 22](../part-7-numerical-computational/22-numerical-methods.md),
- a proof of Euler's formula $e^{j\theta}=\cos\theta+j\sin\theta$.

### Building the Taylor polynomial

We look for a polynomial $p_n(x)=c_0+c_1(x-a)+c_2(x-a)^2+\cdots+c_n(x-a)^n$ with $p_n^{(k)}(a)=f^{(k)}(a)$ for $k=0,1,\ldots,n$. Differentiating $p_n$ $k$ times kills all terms of degree below $k$ and turns $c_k(x-a)^k$ into $k!\,c_k$; the higher terms still contain a factor $(x-a)$ and vanish at $x=a$. So $p_n^{(k)}(a)=k!\,c_k$, and matching gives $c_k=f^{(k)}(a)/k!$.

> [!definition] Taylor polynomial and series
> The **Taylor polynomial of degree $n$** of $f$ about $x=a$ is
> $$
> p_n(x)=\sum_{k=0}^{n}\frac{f^{(k)}(a)}{k!}(x-a)^k=f(a)+f'(a)(x-a)+\frac{f''(a)}{2!}(x-a)^2+\cdots+\frac{f^{(n)}(a)}{n!}(x-a)^n .
> $$
> Letting $n\to\infty$ gives the **Taylor series**. With $a=0$ it is called a **Maclaurin series**.

The degree-1 polynomial $p_1(x)=f(a)+f'(a)(x-a)$ is the tangent line: the linearization of §4.8. Each further term corrects the curvature, then the rate of change of curvature, and so on.

### Taylor's theorem: how big is the error?

> [!theorem] Taylor's theorem with Lagrange remainder
> If $f$ has $n+1$ continuous derivatives on an interval containing $a$ and $x$, then
> $$
> f(x)=p_n(x)+R_n(x),\qquad R_n(x)=\frac{f^{(n+1)}(c)}{(n+1)!}(x-a)^{n+1}
> $$
> for some $c$ between $a$ and $x$. Consequently, if $|f^{(n+1)}|\le M$ between $a$ and $x$,
> $$
> |R_n(x)|\le\frac{M\,|x-a|^{n+1}}{(n+1)!}.
> $$

*Proof sketch.* For $n=0$ this is the mean value theorem: $f(x)=f(a)+f'(c)(x-a)$. For general $n$, fix $x$ and choose the constant $K$ so that $g(t)=f(t)-p_n(t)-K(t-a)^{n+1}$ vanishes at $t=x$. Then $g$ and its first $n$ derivatives vanish at $t=a$, and $g(x)=0$ too. Rolle's theorem gives a zero of $g'$ between $a$ and $x$; applying Rolle repeatedly to $g',g'',\ldots$ gives a point $c$ with $g^{(n+1)}(c)=0$, i.e. $f^{(n+1)}(c)=(n+1)!\,K$. Hence $R_n(x)=K(x-a)^{n+1}$ has the stated form. $\blacksquare$

An equivalent integral form, $R_n(x)=\frac{1}{n!}\int_a^x f^{(n+1)}(t)(x-t)^n\,dt$, follows by repeated integration by parts ([Ch. 5](05-integral-calculus.md)).

The remainder says two things. Near $a$, the error behaves like $|x-a|^{n+1}$: halving the distance cuts the error of a degree-$n$ approximation by $2^{n+1}$. Far from $a$, the factorial in the denominator usually wins in the end, but only if the derivatives do not grow faster.

### The standard series

Computing derivatives at $0$ gives the following Maclaurin series, valid in the stated ranges:

| Function | Maclaurin series | Converges for |
|---|---|---|
| $e^x$ | $1+x+\dfrac{x^2}{2!}+\dfrac{x^3}{3!}+\cdots=\sum\dfrac{x^n}{n!}$ | all $x$ |
| $\sin x$ | $x-\dfrac{x^3}{3!}+\dfrac{x^5}{5!}-\cdots$ | all $x$ |
| $\cos x$ | $1-\dfrac{x^2}{2!}+\dfrac{x^4}{4!}-\cdots$ | all $x$ |
| $\sinh x$, $\cosh x$ | as $\sin$, $\cos$ but all signs $+$ | all $x$ |
| $\dfrac{1}{1-x}$ | $1+x+x^2+x^3+\cdots$ | $\lvert x\rvert<1$ |
| $\ln(1+x)$ | $x-\dfrac{x^2}{2}+\dfrac{x^3}{3}-\cdots$ | $-1<x\le1$ |
| $(1+x)^\alpha$ | $1+\alpha x+\dfrac{\alpha(\alpha-1)}{2!}x^2+\dfrac{\alpha(\alpha-1)(\alpha-2)}{3!}x^3+\cdots$ | $\lvert x\rvert<1$ |
| $\arctan x$ | $x-\dfrac{x^3}{3}+\dfrac{x^5}{5}-\cdots$ | $\lvert x\rvert\le1$ |
| $\tan x$ | $x+\dfrac{x^3}{3}+\dfrac{2x^5}{15}+\cdots$ | $\lvert x\rvert<\pi/2$ |

*Derivations.* For $e^x$, every derivative at $0$ is $1$. For $\sin x$, the derivatives at $0$ cycle $0,1,0,-1$, so only odd powers appear, with alternating signs. For $(1+x)^\alpha$, the $k$-th derivative at $0$ is $\alpha(\alpha-1)\cdots(\alpha-k+1)$; when $\alpha$ is a positive integer the series terminates and is the binomial theorem. The $\ln(1+x)$ and $\arctan x$ series are obtained by integrating $1/(1+x)=1-x+x^2-\cdots$ and $1/(1+x^2)=1-x^2+x^4-\cdots$ term by term.

```desmos-graph
left=-7; right=7; bottom=-3; top=3
---
y=\sin(x)|black
y=x|red|dashed
y=x-\frac{x^3}{6}|orange|dashed
y=x-\frac{x^3}{6}+\frac{x^5}{120}|green|dashed
y=x-\frac{x^3}{6}+\frac{x^5}{120}-\frac{x^7}{5040}|blue|dashed
```
*Maclaurin polynomials of $\sin x$ of degree 1, 3, 5 and 7. Each added term extends the interval over which the polynomial hugs the curve; all of them fail eventually, since a polynomial is unbounded and $\sin x$ is not.*

### Radius of convergence

A power series $\sum c_n(x-a)^n$ converges for $|x-a|<\rho$ and diverges for $|x-a|>\rho$, where the **radius of convergence** $\rho$ can often be found from the ratio test:

$$
\rho=\lim_{n\to\infty}\left|\frac{c_n}{c_{n+1}}\right|\quad(\text{when this limit exists}).
$$

For $e^x$, $c_n=1/n!$ and $c_n/c_{n+1}=n+1\to\infty$: it converges everywhere. For $1/(1-x)$, $c_n=1$ and $\rho=1$. Inside the radius, a power series can be differentiated and integrated term by term, and the result has the same radius.

> [!tip] Engineering Note — why $1/(1+x^2)$ has radius 1
> The function $1/(1+x^2)$ is perfectly smooth for every real $x$, yet its Maclaurin series $1-x^2+x^4-\cdots$ diverges for $|x|>1$. The reason is invisible on the real line: as a function of a complex variable, $1/(1+z^2)$ has poles at $z=\pm j$, a distance $1$ from the origin. **The radius of convergence is the distance to the nearest singularity in the complex plane.** This is why the pole locations of a transfer function control how far a series expansion (a low-frequency approximation, say) can be trusted. The full story is in [Ch. 18](../part-5-transforms/18-complex-analysis.md).

### Euler's formula from series

Substitute $x=j\theta$ into the exponential series and separate even and odd powers, using $j^2=-1$, $j^3=-j$, $j^4=1$:

$$
e^{j\theta}=1+j\theta-\frac{\theta^2}{2!}-j\frac{\theta^3}{3!}+\frac{\theta^4}{4!}+\cdots=\left(1-\frac{\theta^2}{2!}+\frac{\theta^4}{4!}-\cdots\right)+j\left(\theta-\frac{\theta^3}{3!}+\cdots\right)=\cos\theta+j\sin\theta .
$$

The rearrangement is legitimate because the exponential series converges absolutely for every complex argument.

> [!example] Worked Example 4.29 — Computing with a remainder bound
> (a) Estimate $e^{0.1}$ with the degree-3 Maclaurin polynomial and bound the error.
>
> $p_3(0.1)=1+0.1+0.005+0.000166\overline{6}=1.1051667$. The remainder is $R_3=e^c(0.1)^4/4!$ with $0<c<0.1$, and $e^c<e^{0.1}<1.2$, so $|R_3|<1.2\times10^{-4}/24=5\times10^{-6}$.
>
> **Check.** $e^{0.1}=1.1051709$; actual error $4.25\times10^{-6}<5\times10^{-6}$. ✓
>
> (b) Estimate $\sin0.5$ with $p_5$ and bound the error.
>
> $p_5(0.5)=0.5-0.0208333+0.0002604=0.4794271$. Because the $x^6$ coefficient is zero, $p_5=p_6$, so we may use the degree-6 remainder: $|R_6|\le|0.5|^7/7!=1.550\times10^{-6}$ (all derivatives of $\sin$ are bounded by $1$).
>
> **Check.** $\sin0.5=0.4794255$; actual error $1.545\times10^{-6}$, just inside the bound. ✓

### Manipulating series

Finding high derivatives is tedious. In practice series are built from the table by substitution, multiplication, division, differentiation and integration.

- **Substitution.** $e^{-x^2}=1-x^2+\dfrac{x^4}{2!}-\dfrac{x^6}{3!}+\cdots$ (replace $x$ by $-x^2$ in $e^x$). This is how the Gaussian error function is evaluated ([Ch. 19](../part-6-probability-statistics/19-probability.md)).
- **Multiplication.** $e^x\sin x=\left(1+x+\tfrac{x^2}2+\cdots\right)\left(x-\tfrac{x^3}6+\cdots\right)=x+x^2+\tfrac{x^3}{3}+\cdots$.
- **Division.** For $\tan x=\sin x/\cos x$, write $\tan x=b_1x+b_3x^3+b_5x^5+\cdots$ (odd function) and require $(b_1x+b_3x^3+b_5x^5)(1-\tfrac{x^2}2+\tfrac{x^4}{24})=x-\tfrac{x^3}6+\tfrac{x^5}{120}$. Matching powers: $b_1=1$; $b_3-\tfrac12=-\tfrac16\Rightarrow b_3=\tfrac13$; $b_5-\tfrac{b_3}{2}+\tfrac1{24}=\tfrac1{120}\Rightarrow b_5=\tfrac{2}{15}$.
- **Limits.** Series settle indeterminate forms instantly: $\dfrac{x-\sin x}{x^3}=\dfrac{x^3/6-x^5/120+\cdots}{x^3}\to\dfrac16$. (Three applications of L'Hôpital's rule give the same answer more slowly.)

> [!example] Worked Example 4.30 — A Taylor series about a point other than $0$
> Expand $f(x)=\sqrt x$ about $a=4$ to second order and use it to estimate $\sqrt{4.1}$.
>
> $f=x^{1/2}$, $f'=\tfrac12x^{-1/2}$, $f''=-\tfrac14x^{-3/2}$. At $x=4$: $f=2$, $f'=\tfrac14$, $f''=-\tfrac1{32}$. So
> $$
> \sqrt x\approx2+\frac{x-4}{4}-\frac{(x-4)^2}{64}.
> $$
> At $x=4.1$: $2+0.025-0.00015625=2.02484375$.
>
> **Check.** $\sqrt{4.1}=2.02484567$; the error is $1.9\times10^{-6}$, consistent with the next term $f'''(4)(0.1)^3/6=\tfrac{3}{8}\cdot4^{-5/2}\cdot\tfrac{10^{-3}}{6}=1.95\times10^{-6}$. ✓

> [!example] Worked Example 4.31 — Harmonic distortion of an exponential device
> A forward-biased diode (or a bipolar transistor) carries $i=I_Qe^{v/V_T}$, where $v$ is the small signal added to the bias, $V_T=26\,$mV. Let $v=V_m\cos\omega t$. Estimate the second-harmonic distortion $\mathrm{HD}_2$ (amplitude of the $2\omega$ component divided by that of the $\omega$ component) when $V_m=2.6\,$mV.
>
> **Expand** with $x=v/V_T=a\cos\omega t$, $a=V_m/V_T$:
> $$
> \frac{i}{I_Q}=1+a\cos\omega t+\frac{a^2}{2}\cos^2\omega t+\frac{a^3}{6}\cos^3\omega t+\cdots
> $$
> **Use** $\cos^2\theta=\tfrac12(1+\cos2\theta)$ ([Ch. 3](../part-1-foundations/03-trigonometry-and-periodic-functions.md)):
> $$
> \frac{i}{I_Q}\approx\underbrace{1+\frac{a^2}{4}}_{\text{DC shift}}+\underbrace{a\cos\omega t}_{\text{fundamental}}+\underbrace{\frac{a^2}{4}\cos2\omega t}_{\text{2nd harmonic}}+\cdots
> $$
> So $\mathrm{HD}_2\approx\dfrac{a^2/4}{a}=\dfrac a4=\dfrac{V_m}{4V_T}$.
>
> **Numbers.** $a=2.6/26=0.1$, so $\mathrm{HD}_2\approx2.5\,\%$. To keep distortion below $1\,\%$ the input must stay below $4\times0.01\times26=1.04\,$mV peak.
>
> **Check.** The exact Fourier coefficients of $e^{a\cos\omega t}$ are modified Bessel functions, giving $\mathrm{HD}_2=I_2(a)/I_1(a)=0.02499$ at $a=0.1$. ✓ (At $a=0.5$ the estimate $0.125$ is still within $1\,\%$ of the exact $0.1237$.)
>
> The DC shift $a^2/4$ is the mechanism of diode detection and mixing: a nonlinear device converts part of an AC signal into DC.

### Finite-difference formulas from Taylor series

Write the Taylor expansions at $x\pm h$:

$$
f(x\pm h)=f(x)\pm hf'(x)+\frac{h^2}{2}f''(x)\pm\frac{h^3}{6}f'''(x)+\frac{h^4}{24}f^{(4)}(x)\pm\cdots
$$

Subtracting, adding, or using one of them gives the standard numerical derivatives and their errors:

$$
\begin{aligned}
\text{forward:}\quad& \frac{f(x+h)-f(x)}{h}=f'(x)+\frac h2f''(x)+O(h^2),\\
\text{central:}\quad& \frac{f(x+h)-f(x-h)}{2h}=f'(x)+\frac{h^2}{6}f'''(x)+O(h^4),\\
\text{second:}\quad& \frac{f(x+h)-2f(x)+f(x-h)}{h^2}=f''(x)+\frac{h^2}{12}f^{(4)}(x)+O(h^4).
\end{aligned}
$$

For $f=e^x$ at $x=0$ with $h=0.1$, the forward difference is off by $0.0517$ (predicted $h/2\approx0.05$), the central difference by $0.00167$ (predicted $h^2/6=0.00167$) and the second difference by $0.00083$ (predicted $h^2/12=0.00083$). The central difference is far better for the same number of function evaluations, because the even-order error terms cancel by symmetry.

> [!warning] Common Pitfall — using a series outside its radius
> $\ln(1+x)=x-x^2/2+x^3/3-\cdots$ is useless at $x=2$: the terms grow and the partial sums oscillate wildly. To compute $\ln3$, use $\ln3=\ln\frac{1+y}{1-y}$ with $y=\tfrac12$, and the rapidly convergent series $\ln\frac{1+y}{1-y}=2\left(y+\frac{y^3}{3}+\frac{y^5}{5}+\cdots\right)$. Similarly, the binomial approximation $(1+x)^\alpha\approx1+\alpha x$ is only good for $|x|\ll1$, and it is a statement about $|x|$, not $|\alpha x|$.

## 4.8 Linearization

### The tangent-line approximation

> [!definition] Linearization
> The **linearization** of $f$ at $a$ is the degree-1 Taylor polynomial
> $$
> L(x)=f(a)+f'(a)(x-a).
> $$
> For $x$ near $a$, $f(x)\approx L(x)$, with error $f(x)-L(x)=\tfrac12f''(c)(x-a)^2$ for some $c$ between $a$ and $x$.

In terms of increments, with $\Delta x=x-a$ and $\Delta f=f(x)-f(a)$:

$$
\Delta f\approx f'(a)\,\Delta x .
$$

The **differential** notation writes this as $df=f'(x)\,dx$: $dx$ is an independent increment, and $df$ is the corresponding change *along the tangent line*. The true change $\Delta f$ and the differential $df$ differ by a term of order $(\Delta x)^2$.

Linearization is the most important approximation in engineering. Nonlinear devices (diodes, transistors, magnetic cores, sensors) are analysed by choosing an **operating point** (a bias, or quiescent point, Q) and replacing the device by its tangent line there. The slope becomes a small-signal parameter: a conductance, a transconductance, an incremental inductance.

```desmos-graph
left=-0.5; right=9; bottom=-0.5; top=3.5
---
y=\sqrt{x}|blue|x>0
y=2+\frac{x-4}{4}|red|dashed
(4,2)|label:operating point
```
*$\sqrt x$ and its linearization at $x=4$. Near the operating point the line is an excellent substitute; the error grows quadratically with distance, and the line always lies above this concave-down curve.*

> [!example] Worked Example 4.32 — Quick numerical estimates
> (a) $\sqrt{4.1}$: with $f=\sqrt x$, $a=4$, $f'(4)=\tfrac14$: $L(4.1)=2+0.1/4=2.025$. True: $2.024846$ (error $1.5\times10^{-4}$, which matches $\tfrac12|f''(4)|(0.1)^2=1.6\times10^{-4}$).
>
> (b) $\sin31^\circ$: with $f=\sin x$, $a=30^\circ=\pi/6$, $\Delta x=1^\circ=\pi/180=0.017453$ rad: $L=\sin\frac\pi6+\cos\frac\pi6\cdot\frac{\pi}{180}=0.5+0.866\times0.017453=0.51511$. True: $0.51504$. (Converting $\Delta x$ to radians is essential.)
>
> (c) $1.02^{10}$: with $(1+x)^n\approx1+nx$: $1+10(0.02)=1.20$. True: $1.2190$. The error, $\approx\binom{10}{2}(0.02)^2=0.018$, is the next binomial term.
>
> (d) Small angles: $\sin\theta\approx\theta$ at $\theta=10^\circ=0.17453$ rad gives $0.17453$ against the true $0.17365$, an error of $0.5\,\%$. The relative error of $\sin\theta\approx\theta$ is about $\theta^2/6$.

### Small-signal models

> [!example] Worked Example 4.33 — MOSFET transconductance
> A MOSFET in saturation obeys the square law $i_D=\tfrac{k}{2}(v_{GS}-V_t)^2$ with $k=2\,$mA/V² and $V_t=0.5\,$V. It is biased at $V_{GS}=1.0\,$V. Find $I_D$, the transconductance $g_m=di_D/dv_{GS}$ at the Q-point, and the error of the linear model for a $10\,$mV signal.
>
> **Bias.** Overdrive $V_{ov}=V_{GS}-V_t=0.5\,$V, so $I_D=\tfrac{2}{2}(0.5)^2=0.25\,$mA.
>
> **Linearize.** $g_m=k(V_{GS}-V_t)=2\times0.5=1\,$mA/V. Equivalently $g_m=2I_D/V_{ov}=2(0.25)/0.5=1\,$mA/V.
>
> **Small-signal model.** $i_D\approx I_D+g_mv_{gs}$, where $v_{gs}=v_{GS}-V_{GS}$ is the signal. For $v_{gs}=10\,$mV, the model predicts $\Delta i_D=10\,\mu$A.
>
> **Exact.** $\Delta i_D=\tfrac k2\big[(0.51)^2-(0.5)^2\big]=10^{-3}(0.2601-0.25)=10.1\,\mu$A. The relative error is $v_{gs}/(2V_{ov})=1\,\%$. Because the square law has $f''$ constant, the error is exactly the quadratic term $\tfrac k2v_{gs}^2=0.1\,\mu$A.

> [!example] Worked Example 4.34 — Linearizing in two variables
> The magnitude of an impedance $Z=R+jX$ is $|Z|=\sqrt{R^2+X^2}$. Linearize it about $R=3\,\Omega$, $X=4\,\Omega$ and estimate $|Z|$ at $R=3.1$, $X=3.9$.
>
> The two-variable linearization uses both partial derivatives (§4.5):
> $$
> L(R,X)=|Z|_0+\frac{\partial|Z|}{\partial R}\Delta R+\frac{\partial|Z|}{\partial X}\Delta X,
> \qquad \frac{\partial|Z|}{\partial R}=\frac{R}{|Z|},\ \frac{\partial|Z|}{\partial X}=\frac{X}{|Z|}.
> $$
> At $(3,4)$: $|Z|_0=5$, partials $0.6$ and $0.8$. So $L=5+0.6(0.1)+0.8(-0.1)=4.98\,\Omega$.
>
> **Check.** $\sqrt{3.1^2+3.9^2}=\sqrt{24.82}=4.98197\,\Omega$; error $0.002\,\Omega$. ✓ Geometrically, $L(R,X)$ is the tangent plane to the cone $z=\sqrt{R^2+X^2}$ ([Ch. 6](06-multivariable-calculus.md)).

### Linearizing a differential equation about an equilibrium

Consider a first-order system $\dot x=f(x)$. An **equilibrium** $x^*$ satisfies $f(x^*)=0$. Write $x=x^*+\delta$ with $\delta$ small. Since $\dot x^*=0$,

$$
\dot\delta=f(x^*+\delta)\approx f(x^*)+f'(x^*)\,\delta=f'(x^*)\,\delta .
$$

This linear equation has solution $\delta(t)=\delta(0)e^{f'(x^*)t}$. Hence:

- if $f'(x^*)<0$, small disturbances decay: the equilibrium is **stable**, with time constant $\tau=-1/f'(x^*)$;
- if $f'(x^*)>0$, they grow: the equilibrium is **unstable**;
- if $f'(x^*)=0$, linearization is inconclusive.

For systems of equations the derivative becomes a Jacobian matrix and its eigenvalues decide stability ([Ch. 10](../part-3-linear-algebra/10-eigenvalues-and-eigenvectors.md)).

> [!example] Worked Example 4.35 — A capacitor charged by a current source into a diode
> A constant current $I_0=1\,$mA charges a capacitor $C=1\,$nF that is shunted by a diode ($I_s=10^{-14}\,$A, $V_T=25.85\,$mV). The capacitor voltage obeys
> $$
> C\frac{dv}{dt}=I_0-I_s\big(e^{v/V_T}-1\big)\equiv C\,f(v).
> $$
> Find the equilibrium, decide its stability and find the time constant of small deviations.
>
> **Equilibrium.** $I_s(e^{v^*/V_T}-1)=I_0\Rightarrow v^*=V_T\ln(1+I_0/I_s)=0.02585\ln(10^{11}+1)=0.6547\,$V. (All of $I_0$ flows through the diode.)
>
> **Linearize.** $f'(v)=-\dfrac{I_s}{CV_T}e^{v/V_T}$, so $f'(v^*)=-\dfrac{I_0+I_s}{CV_T}<0$: stable.
>
> **Time constant.** $\tau=-1/f'(v^*)=\dfrac{CV_T}{I_0+I_s}\approx C\,r_d$ with $r_d=V_T/I_0=25.85\,\Omega$. Numerically $\tau=10^{-9}\times25.85=25.9\,$ns.
>
> The linearized circuit is the capacitor in parallel with the diode's small-signal resistance: an RC circuit. Linearization has turned a nonlinear ODE into the most familiar circuit of all.

### Newton's method

Linearization also solves equations. To solve $f(x)=0$ starting from a guess $x_n$, replace $f$ by its tangent line at $x_n$ and solve the *linear* equation $f(x_n)+f'(x_n)(x-x_n)=0$:

$$
x_{n+1}=x_n-\frac{f(x_n)}{f'(x_n)} .
$$

> [!theorem] Quadratic convergence of Newton's method
> If $f$ is $C^2$ near a simple root $r$ ($f(r)=0$, $f'(r)\neq0$) and $x_0$ is close enough to $r$, the errors $e_n=x_n-r$ satisfy
> $$
> e_{n+1}\approx\frac{f''(r)}{2f'(r)}\,e_n^2 .
> $$
> The number of correct digits roughly doubles with each step.

*Proof sketch.* Taylor's theorem about $x_n$: $0=f(r)=f(x_n)+f'(x_n)(r-x_n)+\tfrac12f''(c)(r-x_n)^2$. Divide by $f'(x_n)$ and rearrange: $x_{n+1}-r=\dfrac{f''(c)}{2f'(x_n)}(x_n-r)^2$.

> [!example] Worked Example 4.36 — Newton's method for a diode bias point
> A $5\,$V source drives a $1\,$kΩ resistor in series with a diode ($I_s=10^{-14}\,$A, $V_T=25.852\,$mV). Find the diode voltage $V$.
>
> KCL at the diode: $f(V)=\dfrac{5-V}{1000}-I_s\big(e^{V/V_T}-1\big)=0$, with $f'(V)=-\dfrac{1}{1000}-\dfrac{I_s}{V_T}e^{V/V_T}$.
>
> Starting from $V_0=0.7\,$V:
>
> | $n$ | $V_n$ (V) | error $V_n-V^*$ |
> |---|---|---|
> | 0 | 0.700000 | $7.5\times10^{-3}$ |
> | 1 | 0.693518 | $9.7\times10^{-4}$ |
> | 2 | 0.692562 | $1.8\times10^{-5}$ |
> | 3 | 0.6925436394 | $6.3\times10^{-9}$ |
> | 4 | 0.6925436332 | $<10^{-15}$ |
>
> The diode voltage is $V^*=0.69254\,$V and the current $(5-0.69254)/1000=4.3075\,$mA. The error is roughly squared at each step, as the theorem predicts.

> [!warning] Common Pitfall — Newton's method far from the root
> Started at $V_0=0\,$V, the iteration above jumps to $V_1\approx5.0\,$V (at $0\,$V the diode is nearly an open circuit, so the tangent line only sees the resistor) and then creeps back down by about one $V_T$ per step: it needs roughly 170 iterations, and along the way $e^{V/V_T}\approx e^{193}\approx10^{84}$. With a slightly larger source voltage the exponential overflows in floating point. Circuit simulators therefore **limit** the voltage step per iteration on junctions. Newton's method is excellent locally and unreliable globally; combine it with bracketing (§4.1) or step limiting.

> [!tip] Engineering Note — the limits of small-signal analysis
> A linear model is valid only while the neglected quadratic term is small compared with the linear one: for the exponential device, $|v|\ll2V_T\approx50\,$mV; for the square-law MOSFET, $|v_{gs}|\ll2V_{ov}$. Example 4.31 quantified exactly what "small" costs in distortion. Large-signal behaviour (clipping, slewing, switching) requires the full nonlinear model, usually solved numerically.

## 4.9 Rates of change

### Average and instantaneous rates

For a quantity $y=f(t)$, the **average rate of change** over $[t_1,t_2]$ is $\dfrac{f(t_2)-f(t_1)}{t_2-t_1}$, the slope of a secant. The **instantaneous rate of change** at $t$ is $f'(t)$, the slope of the tangent. The mean value theorem (§4.2) links them: somewhere in the interval the instantaneous rate equals the average.

Almost every constitutive law in circuit theory is a rate law:

| Law | Rate | Units |
|---|---|---|
| Current | $i=\dfrac{dq}{dt}$ | C/s = A |
| Capacitor | $i=C\dfrac{dv}{dt}$ (more generally $i=\dfrac{d(Cv)}{dt}$) | F·V/s = A |
| Inductor | $v=L\dfrac{di}{dt}$ (more generally $v=\dfrac{d\lambda}{dt}$, $\lambda$ = flux linkage) | H·A/s = V |
| Faraday's law | $\mathcal{E}=-\dfrac{d\Phi}{dt}$ | Wb/s = V |
| Power | $p=\dfrac{dW}{dt}$ | J/s = W |
| Slew rate | $\mathrm{SR}=\max\left\lvert\dfrac{dv_o}{dt}\right\rvert$ | V/s (usually V/µs) |

> [!example] Worked Example 4.37 — Average versus instantaneous current
> The charge through a conductor is $q(t)=3t^2-t^3\,$mC for $0\le t\le2\,$s. Find the average current over $[0,2]$, the instantaneous current at $t=1\,$s, and the times at which the instantaneous current equals the average.
>
> **Average.** $\bar i=\dfrac{q(2)-q(0)}{2}=\dfrac{12-8}{2}=2\,$mA.
>
> **Instantaneous.** $i(t)=q'(t)=6t-3t^2\,$mA, so $i(1)=3\,$mA (the maximum, since $i'(t)=6-6t=0$ there).
>
> **MVT points.** $6c-3c^2=2\Rightarrow3c^2-6c+2=0\Rightarrow c=1\pm\dfrac{1}{\sqrt3}=0.423\,$s or $1.577\,$s. Both lie in $(0,2)$, as the MVT guarantees at least one must.

```desmos-graph
left=-0.2; right=2.4; bottom=-0.5; top=5
---
y=3x^2-x^3|blue|0<x<2
y=2x|grey|dashed|0<x<2
y=2(x-0.4226)+3(0.4226)^2-0.4226^3|red|0<x<1.2
y=2(x-1.5774)+3(1.5774)^2-1.5774^3|red|1<x<2.2
```
*The charge $q(t)=3t^2-t^3$ (blue), the secant over $[0,2]$ with slope $2$ (the average current, dashed) and the two tangent lines with the same slope (red) at $t=1\mp1/\sqrt3$.*

### Relative (logarithmic) rates

The **relative rate of change** is $\dfrac{\dot y}{y}=\dfrac{d}{dt}\ln y$, the fraction by which $y$ changes per unit time. For an exponential $y=Y_0e^{-t/\tau}$ it is constant: $\dot y/y=-1/\tau$. That is the defining property of exponential decay: the *relative* rate never changes, so the capacitor in an RC circuit loses the same *fraction* of its remaining voltage in every equal time interval. With $\tau=2\,$ms, the instantaneous relative rate is $-500\,\text{s}^{-1}$, i.e. $-0.05\,\%$ per µs. Note that the voltage falls by $1-e^{-1/2}=39\,\%$, not $50\,\%$, in the first millisecond: a relative rate is instantaneous, not a percentage per finite interval.

### Related rates

In a **related-rates** problem, several quantities vary with time and are linked by an equation. Differentiating the equation with respect to $t$ (with the chain rule, as in implicit differentiation) links their rates. The procedure:

1. Name the variables and write the equation that connects them, valid at **all** times.
2. Differentiate both sides with respect to $t$.
3. Only **then** substitute the values at the instant of interest, and solve for the unknown rate.

Substituting before differentiating is the classic mistake: it freezes quantities that are in fact changing.

> [!example] Worked Example 4.38 — A MEMS capacitive sensor
> A parallel-plate capacitor has $C=\varepsilon_0A/d$ with $A=1\,\text{cm}^2$, held at a constant $V=5\,$V. At the instant when $d=10\,\mu$m, the moving plate recedes at $\dot d=1\,$mm/s. Find the current.
>
> **Step 1.** $q=Cv$ and $C=\varepsilon_0A/d$.
>
> **Step 2.** Product rule: $i=\dot q=\dot Cv+C\dot v$. Here $\dot v=0$. Chain rule: $\dot C=-\dfrac{\varepsilon_0A}{d^2}\dot d=-\dfrac{C}{d}\dot d$.
>
> **Step 3.** $C=8.854\times10^{-12}\times10^{-4}/10^{-5}=88.54\,$pF. $\dot C=-88.54\times10^{-12}\times(10^{-3}/10^{-5})=-8.854\,$nF/s. So
> $$
> i=V\dot C=5\times(-8.854\times10^{-9})=-44.3\,\text{nA}.
> $$
> The negative sign means charge flows *out* of the capacitor as the plates separate. This current is the output signal of a capacitive microphone or accelerometer.
>
> **Units.** $\text{V}\cdot\text{F/s}=\text{C/s}=\text{A}$. ✓

> [!example] Worked Example 4.39 — EMF in an expanding loop
> A circular conducting loop lies in a uniform magnetic field $B=0.5\,$T perpendicular to its plane. Its radius increases at $\dot r=2\,$cm/s. Find the induced EMF when $r=10\,$cm.
>
> **Equation.** $\Phi=B\pi r^2$, valid at all times.
>
> **Differentiate.** $\dot\Phi=2\pi Br\dot r$ ($B$ constant), so $|\mathcal{E}|=2\pi Br\dot r$.
>
> **Substitute.** $|\mathcal{E}|=2\pi(0.5)(0.1)(0.02)=6.28\,$mV.
>
> **Interpretation.** $2\pi r\dot r$ is the rate at which area is added (the circumference times the outward speed), and $B$ times that is the rate at which flux is swept: the motional EMF.

> [!example] Worked Example 4.40 — Slew-rate limit of an amplifier
> An op-amp has slew rate $\mathrm{SR}=0.5\,\text{V}/\mu\text{s}$. What is the highest frequency at which it can output an undistorted $10\,$V-peak sine wave?
>
> For $v_o=V_m\sin(2\pi ft)$, $\dot v_o=2\pi fV_m\cos(2\pi ft)$, whose maximum magnitude is $2\pi fV_m$ (at the zero crossings, where the sine is steepest). Require $2\pi fV_m\le\mathrm{SR}$:
> $$
> f_{\max}=\frac{\mathrm{SR}}{2\pi V_m}=\frac{0.5\times10^6}{2\pi\times10}=7.96\ \text{kHz}.
> $$
> This is the **full-power bandwidth**. Above it the output becomes a triangle wave, even if the amplifier's small-signal bandwidth is far higher. Note that the limit depends on amplitude: at $1\,$V peak, $f_{\max}=79.6\,$kHz.

> [!example] Worked Example 4.41 — Power in reactive elements
> Show that the rate of change of stored energy in an inductor equals the power delivered to it.
>
> $W=\tfrac12Li^2$. Chain rule: $\dfrac{dW}{dt}=\tfrac12L\cdot2i\dfrac{di}{dt}=i\left(L\dfrac{di}{dt}\right)=iv=p$. ✓ Likewise for a capacitor, $W=\tfrac12Cv^2$ gives $\dot W=Cv\dot v=vi$.
>
> For $i=I_m\sin\omega t$: $p=\tfrac12\omega LI_m^2\sin2\omega t$, which averages to zero. An ideal inductor stores and returns energy; it does not dissipate it.

> [!warning] Common Pitfall — signs, units and constants that are not constant
> Decide a sign convention first (passive sign convention: current entering the $+$ terminal) and keep it; a negative rate is information, not an error. Convert units before substituting ($1\,$mm/s $=10^{-3}\,$m/s, $\mu$s in the denominator of a slew rate means $\times10^{6}$). And remember that "$i=C\,dv/dt$" assumes $C$ is constant: for a varying capacitance, as in Example 4.38, the full rule is $i=d(Cv)/dt$.

> [!tip] Engineering Note — rates as specifications
> Many specifications are bounds on derivatives: slew rate (V/µs), $dv/dt$ immunity of thyristors and gate drivers (V/ns), $di/dt$ limits on power semiconductors (A/µs), ramp rates of power-supply sequencing, and the rate of change of frequency (RoCoF, Hz/s) in power grids. The inductor law also explains why a fast current edge through a small parasitic inductance causes trouble: $1\,$A switched in $1\,$ns through $10\,$nH of trace inductance produces $L\,di/dt=10\,$V of ground bounce.

## 4.10 Engineering optimization

### Motivation

Design is optimization: the load that extracts the most power, the switching frequency that minimizes total loss, the cable geometry with least attenuation, the damping that gives the flattest response. When a performance measure depends smoothly on one design variable, calculus finds the best value exactly. Multivariable and constrained problems (Lagrange multipliers) are in [Ch. 6](06-multivariable-calculus.md); numerical optimization algorithms are in [Ch. 23](../part-7-numerical-computational/23-optimization.md).

### Critical points

> [!definition] Extrema and critical points
> $f$ has a **local maximum** at $c$ if $f(c)\ge f(x)$ for all $x$ near $c$, and an **absolute (global) maximum** on a set $D$ if $f(c)\ge f(x)$ for all $x\in D$; similarly for minima. A **critical point** of $f$ is an interior point $c$ of its domain where $f'(c)=0$ or $f'(c)$ does not exist.

> [!theorem] Fermat's theorem
> If $f$ has a local extremum at an interior point $c$ and $f'(c)$ exists, then $f'(c)=0$.

*Proof.* Suppose $c$ is a local maximum. For small $h>0$, $f(c+h)-f(c)\le0$, so $\frac{f(c+h)-f(c)}{h}\le0$ and the right-hand limit gives $f'(c)\le0$. For small $h<0$ the quotient is $\ge0$, giving $f'(c)\ge0$. Hence $f'(c)=0$. $\blacksquare$

The converse fails: $f=x^3$ has $f'(0)=0$ but no extremum at $0$. Critical points are *candidates*, which must then be classified.

> [!theorem] Classifying critical points
> **First-derivative test.** If $f'$ changes sign from $+$ to $-$ at $c$, $f$ has a local maximum at $c$; from $-$ to $+$, a local minimum; no sign change, no extremum.
>
> **Second-derivative test.** If $f'(c)=0$ and $f''(c)<0$, $c$ is a local maximum; if $f''(c)>0$, a local minimum; if $f''(c)=0$, the test is inconclusive.

*Why the second-derivative test works.* By Taylor's theorem, $f(c+h)=f(c)+0\cdot h+\tfrac12f''(c)h^2+o(h^2)$. For small $h$ the quadratic term dominates, and its sign is that of $f''(c)$, whichever sign $h$ has.

### The closed-interval method

When the design variable lies in a closed interval $[a,b]$ (a component range, a physical limit), the EVT guarantees an absolute maximum and minimum, and they occur either at critical points or at the endpoints. So:

1. Find the critical points in $(a,b)$.
2. Evaluate $f$ at them and at $a$ and $b$.
3. The largest value is the absolute maximum, the smallest the absolute minimum.

For example, $f(x)=x^3-3x$ on $[-2,3]$ has critical points $x=\pm1$ (from $3x^2-3=0$). The values are $f(-2)=-2$, $f(-1)=2$, $f(1)=-2$, $f(3)=18$. The absolute maximum $18$ is at the endpoint $x=3$, not at the local maximum $x=-1$; the absolute minimum $-2$ is attained twice.

### A strategy for applied problems

1. Identify the design variable and the objective, with units; write the objective as a function of one variable (use constraints to eliminate others).
2. State the feasible range (resistances are positive, frequencies lie in a band).
3. Differentiate, solve $f'=0$, and classify.
4. Compare with endpoints and check that the answer is physically sensible.
5. Ask how sensitive the optimum is (§4.11): a flat optimum is forgiving.

> [!example] Worked Example 4.42 — Maximum power transfer
> Choose $R_L>0$ to maximize $P=\dfrac{V^2R_L}{(R_s+R_L)^2}$.
>
> From Example 4.8, $P'(R_L)=\dfrac{V^2(R_s-R_L)}{(R_s+R_L)^3}$. The only critical point in $(0,\infty)$ is $R_L=R_s$. The derivative changes from $+$ to $-$ there, so it is a maximum; since $P\to0$ as $R_L\to0$ and as $R_L\to\infty$, it is the global maximum.
> $$
> P_{\max}=\frac{V^2R_s}{(2R_s)^2}=\frac{V^2}{4R_s}.
> $$
> **Numbers.** $V=10\,$V, $R_s=50\,\Omega$: $P_{\max}=0.5\,$W at $R_L=50\,\Omega$.
>
> **Efficiency.** At the optimum, the source resistance dissipates as much as the load: efficiency is $50\,\%$. Maximum power transfer is the right objective for weak signals (antennas, sensors, RF stages), never for power distribution, where $R_L\gg R_s$ is used to maximize efficiency.
>
> **AC version.** For a source impedance $Z_s=R_s+jX_s$ and load $Z_L=R_L+jX_L$, $P=\tfrac12|V|^2R_L/\big[(R_s+R_L)^2+(X_s+X_L)^2\big]$. Setting $\partial P/\partial X_L=0$ gives $X_L=-X_s$; then $\partial P/\partial R_L=0$ gives $R_L=R_s$. So $Z_L=Z_s^*$: **conjugate matching**.

```desmos-graph
left=0; right=300; bottom=0; top=0.6
---
y=\frac{100x}{(50+x)^2}|blue|x>0
x=50|red|dashed
(50,0.5)|label:max 0.5 W at R_L = R_s
```
*Load power versus load resistance for $V=10\,$V, $R_s=50\,\Omega$. The peak is broad: the curve is flat near the optimum because $P'=0$ there, which is the subject of Example 4.49.*

> [!example] Worked Example 4.43 — Optimal switching frequency
> In a switching converter, losses that scale with ripple (magnetics, capacitor ESR) fall as $1/f$ while switching losses rise as $f$. A model is $P(f)=\dfrac af+bf$ with $a=2\times10^5\,$W·Hz and $b=2\times10^{-6}\,$W/Hz. Find the frequency of minimum loss.
>
> $P'(f)=-\dfrac{a}{f^2}+b=0\Rightarrow f^*=\sqrt{a/b}$. $P''(f)=2a/f^3>0$, so it is a minimum (and the only critical point for $f>0$, with $P\to\infty$ at both ends, so it is global).
> $$
> f^*=\sqrt{\frac{2\times10^5}{2\times10^{-6}}}=\sqrt{10^{11}}=316\ \text{kHz},\qquad P(f^*)=2\sqrt{ab}=2\sqrt{0.4}=1.26\ \text{W}.
> $$
> At the optimum the two loss components are equal ($a/f^*=bf^*=\sqrt{ab}$). This **balance principle** holds whenever a cost is the sum of a term growing like $x$ and one falling like $1/x$: copper versus core loss in a transformer, conduction versus switching loss, wire size versus cost.

> [!example] Worked Example 4.44 — Resonant peak of a second-order low-pass filter
> A second-order low-pass filter has
> $$
> |H(u)|=\frac{1}{\sqrt{(1-u^2)^2+(2\zeta u)^2}},\qquad u=\omega/\omega_0,
> $$
> with damping ratio $\zeta>0$. Find the frequency and height of the resonant peak.
>
> **Simplify first.** $|H|$ is largest where $D(w)=(1-w)^2+4\zeta^2w$ is smallest, with $w=u^2\ge0$. (Maximizing a positive function is equivalent to minimizing its reciprocal; squaring is monotonic for positive quantities, so working with $D$ avoids square roots.)
>
> **Differentiate.** $D'(w)=-2(1-w)+4\zeta^2=0\Rightarrow w^*=1-2\zeta^2$. $D''=2>0$: minimum.
>
> **Feasibility.** $w^*\ge0$ requires $\zeta\le1/\sqrt2$. For $\zeta>1/\sqrt2$, $D$ is increasing on $w\ge0$, the minimum is at the endpoint $w=0$, and there is no peak: the response falls monotonically. $\zeta=1/\sqrt2$ is the **Butterworth** (maximally flat) boundary.
>
> **Peak.** For $\zeta<1/\sqrt2$:
> $$
> u_{\text{peak}}=\sqrt{1-2\zeta^2},\qquad D(w^*)=4\zeta^2-4\zeta^4=4\zeta^2(1-\zeta^2),\qquad |H|_{\max}=\frac{1}{2\zeta\sqrt{1-\zeta^2}}.
> $$
> **Numbers.** $\zeta=0.2$: $u_{\text{peak}}=0.959$ and $|H|_{\max}=2.55$ ($+8.1\,$dB). $\zeta=0.5$: $u_{\text{peak}}=0.707$, $|H|_{\max}=1.155$ ($+1.25\,$dB).

```desmos-graph
left=0; right=2.5; bottom=0; top=3
---
y=\frac{1}{\sqrt{(1-x^2)^2+(0.4x)^2}}|blue
y=\frac{1}{\sqrt{(1-x^2)^2+(1.0x)^2}}|green
y=\frac{1}{\sqrt{(1-x^2)^2+(1.414x)^2}}|red
(0.959,2.55)|label:zeta 0.2 peak
```
*Second-order low-pass magnitude for $\zeta=0.2$ (blue), $0.5$ (green) and $1/\sqrt2$ (red, Butterworth). The peak moves below $\omega_0$ and disappears at $\zeta=1/\sqrt2$, as found in Example 4.44.*

> [!example] Worked Example 4.45 — Why coaxial cable is $75\,\Omega$
> For a coaxial line with inner radius $a$ and fixed outer radius $b$ (set by the cable size), the conductor loss per unit length is proportional to
> $$
> \alpha(a)\propto\frac{1/a+1/b}{\ln(b/a)} .
> $$
> Find the ratio $x=b/a$ that minimizes the loss.
>
> **One variable.** With $1/a=x/b$, $\alpha\propto\dfrac1b\cdot\dfrac{x+1}{\ln x}$; minimize $g(x)=\dfrac{x+1}{\ln x}$ for $x>1$.
>
> **Differentiate.** Quotient rule:
> $$
> g'(x)=\frac{\ln x-(x+1)/x}{(\ln x)^2}=0\quad\Longrightarrow\quad\ln x=1+\frac1x .
> $$
> This transcendental equation has no closed-form solution. Newton's method (§4.8) on $h(x)=\ln x-1-1/x$, $h'(x)=1/x+1/x^2$, from $x_0=3$: $x_1=3.5281$, $x_2=3.5904$, $x_3=3.59112$. So $x^*=3.5911$.
>
> **Check it is a minimum.** $g\to\infty$ as $x\to1^+$ (denominator $\to0$) and as $x\to\infty$ (since $x/\ln x\to\infty$), and $x^*$ is the only critical point, so it is the global minimum.
>
> **Consequence.** The characteristic impedance of an air-filled coax is $Z_0=60\ln(b/a)\,\Omega=60\ln3.5911=76.7\,\Omega$. That is why low-loss video and RF distribution cable is $75\,\Omega$. (The $50\,\Omega$ standard is a compromise between this and the $30\,\Omega$ ratio that maximizes power handling.)

> [!warning] Common Pitfall — forgetting endpoints and the domain
> A critical point is only a candidate. The global optimum may be at an endpoint (as for $x^3-3x$ on $[-2,3]$ above), may not exist on an open interval (minimize $1/x$ on $(0,\infty)$), or the critical point may be infeasible (the resonant peak of Example 4.44 for $\zeta>1/\sqrt2$). Points where $f'$ fails to exist are critical points too: the minimum of $|x-2|$ is at the corner $x=2$. Finally, $f''(c)=0$ does not mean "no extremum": $x^4$ has a minimum at $0$. Use the first-derivative test when the second-derivative test is inconclusive.

> [!tip] Engineering Note — optimizing a simpler function
> Because $\sqrt{\cdot}$, $\ln$, squaring of positive quantities and $x\mapsto1/x$ (on positive values, reversing the order) are monotonic, the optimizer of $f$ is also the optimizer of $f^2$, $\ln f$ or $1/f$. Maximize $|H|^2$ instead of $|H|$, $\ln L$ instead of a likelihood $L$ ([Ch. 20](../part-6-probability-statistics/20-statistical-analysis.md)), and gain in dB instead of linear gain. The location of the optimum is unchanged, and the algebra is usually much lighter.

## 4.11 Sensitivity analysis

### Motivation

Real components are not their nominal values. Resistors come in $1\,\%$ or $5\,\%$ tolerances, capacitors in $10\,\%$ or $20\,\%$; values drift with temperature and age; an op-amp's open-loop gain varies by a factor of several between units. A designer must answer: *how much does the performance change when a component changes?* The answer is a derivative, usually normalized so that percentage changes map to percentage changes. Sensitivity analysis tells you which components need tight tolerances, which can be cheap, and whether a design is robust.

### Absolute and relative sensitivity

> [!definition] Sensitivity
> Let a performance measure $y$ depend on a parameter $x$ (other parameters fixed). The **absolute sensitivity** is $\dfrac{\partial y}{\partial x}$ (units of $y$ per unit of $x$). The **relative (normalized, Bode) sensitivity** is
> $$
> S_x^y=\frac{x}{y}\frac{\partial y}{\partial x}=\frac{\partial y/y}{\partial x/x}=\frac{\partial\ln y}{\partial\ln x}.
> $$
> It is dimensionless: for small changes, a $1\,\%$ change in $x$ produces approximately an $S_x^y\,\%$ change in $y$.

The first-order change (from the linearization of §4.8, in relative form) is

$$
\frac{\Delta y}{y}\approx S_x^y\,\frac{\Delta x}{x}.
$$

Because $S_x^y$ is a logarithmic derivative (§4.3), it obeys simple algebraic rules that make most sensitivities computable at sight.

> [!theorem] Properties of relative sensitivity
> For functions $u(x)$, $v(x)$ and constants $c$, $n$:
> $$
> \begin{aligned}
> S_x^{cu}&=S_x^{u}, & S_x^{uv}&=S_x^u+S_x^v, & S_x^{u/v}&=S_x^u-S_x^v,\\
> S_x^{u^n}&=nS_x^u, & S_x^{x^n}&=n, & S_x^{1/u}&=-S_x^u,\\
> S_x^{y}&=S_u^{y}\,S_x^{u}\ \text{(chain rule)}, & S_x^{u+v}&=\frac{u\,S_x^u+v\,S_x^v}{u+v}.
> \end{aligned}
> $$

*Proofs.* Each follows from $S_x^y=d\ln y/d\ln x$: $\ln(uv)=\ln u+\ln v$ gives the product rule; $\ln u^n=n\ln u$ gives the power rule; and the chain rule is $\dfrac{d\ln y}{d\ln x}=\dfrac{d\ln y}{d\ln u}\dfrac{d\ln u}{d\ln x}$. For the sum, $S_x^{u+v}=\dfrac{x(u'+v')}{u+v}=\dfrac{u\cdot xu'/u+v\cdot xv'/v}{u+v}$.

So for any **monomial formula** $y=c\,x_1^{n_1}x_2^{n_2}\cdots$, the relative sensitivity to $x_k$ is just the exponent $n_k$. For example, for $f_c=\dfrac{1}{2\pi RC}=\dfrac{1}{2\pi}R^{-1}C^{-1}$: $S_R^{f_c}=S_C^{f_c}=-1$.

### Several parameters: worst case and statistical tolerance

If $y$ depends on $x_1,\ldots,x_n$, the total differential of §4.5 gives

$$
\frac{\Delta y}{y}\approx\sum_{k=1}^{n}S_{x_k}^{y}\,\frac{\Delta x_k}{x_k}.
$$

Two ways of combining tolerances follow.

- **Worst case.** If each relative error is bounded, $|\Delta x_k/x_k|\le t_k$, then
$$
\left|\frac{\Delta y}{y}\right|\lesssim\sum_k\big|S_{x_k}^y\big|\,t_k .
$$
This is guaranteed (to first order) but pessimistic: it assumes every component sits at the worst extreme of its tolerance simultaneously.

- **Statistical (root-sum-square, RSS).** If the parameter errors are independent random variables with relative standard deviations $\sigma_k$, then the relative standard deviation of $y$ is approximately
$$
\sigma_y\approx\sqrt{\sum_k\big(S_{x_k}^y\big)^2\sigma_k^2}.
$$
This is the first-order **propagation of uncertainty** ([Ch. 20](../part-6-probability-statistics/20-statistical-analysis.md)). It describes the spread of a production batch, and is smaller than the worst case.

Both rest on linearization, so they are reliable only when the relative changes are small enough that second-order terms are negligible.

> [!example] Worked Example 4.46 — Voltage divider
> $V_o=V_i\dfrac{R_2}{R_1+R_2}$ with $R_1=R_2=10\,$kΩ, both $\pm1\,\%$. Find the sensitivities and the tolerance of $V_o$.
>
> **Sensitivities.** Using the quotient and sum rules, or differentiating directly:
> $$
> S_{R_2}^{V_o}=1-\frac{R_2}{R_1+R_2}=\frac{R_1}{R_1+R_2},\qquad S_{R_1}^{V_o}=-\frac{R_1}{R_1+R_2},\qquad S_{V_i}^{V_o}=1.
> $$
> For equal resistors, $S_{R_2}=+0.5$ and $S_{R_1}=-0.5$.
>
> **Tolerance** (taking $V_i$ as exact). Worst case: $0.5(1\,\%)+0.5(1\,\%)=1\,\%$. RSS: $\sqrt{0.5^2+0.5^2}\,(1\,\%)=0.71\,\%$.
>
> **Exact worst case.** $R_1=9.9\,$kΩ, $R_2=10.1\,$kΩ gives $V_o/V_i=0.505$, i.e. $+1.00\,\%$. ✓
>
> **Insight.** $S_{R_1}+S_{R_2}=0$: the ratio is unchanged if both resistors drift by the same fraction (for example, with a common temperature coefficient). This is why precision dividers use matched resistor networks: tracking matters more than absolute accuracy.

> [!example] Worked Example 4.47 — Cutoff frequency of an RC filter
> $f_c=1/(2\pi RC)$ with $R=10\,$kΩ $\pm5\,\%$ and $C=10\,$nF $\pm10\,\%$ (nominal $f_c=1.592\,$kHz). Estimate the spread of $f_c$.
>
> $S_R^{f_c}=S_C^{f_c}=-1$. Worst case: $5\,\%+10\,\%=15\,\%$. RSS: $\sqrt{5^2+10^2}\,\%=11.2\,\%$.
>
> **Exact extremes.** $R$ and $C$ both high: $f_c/f_{c0}=1/(1.05\times1.10)=0.866$ ($-13.4\,\%$). Both low: $1/(0.95\times0.90)=1.170$ ($+17.0\,\%$). The linear estimate $\pm15\,\%$ is the average of these; the asymmetry is the second-order effect of the reciprocal. For tolerances this large, check the extremes exactly.
>
> **Design consequence.** The capacitor dominates the error. Buying a $2\,\%$ capacitor (RSS $5.4\,\%$) helps far more than buying a $1\,\%$ resistor (RSS $10.05\,\%$).

> [!example] Worked Example 4.48 — Series RLC: resonant frequency and Q
> A series RLC circuit has $\omega_0=\dfrac{1}{\sqrt{LC}}$ and $Q=\dfrac1R\sqrt{\dfrac LC}$. With $L=10\,$mH $\pm5\,\%$, $C=100\,$nF $\pm10\,\%$, $R=10\,\Omega\pm1\,\%$, find the sensitivities and tolerances.
>
> **Nominal values.** $\omega_0=1/\sqrt{10^{-9}}=31\,623\,$rad/s ($f_0=5.03\,$kHz); $Q=0.1\sqrt{10^5}=31.6$.
>
> **Sensitivities (read off the exponents).** $\omega_0=L^{-1/2}C^{-1/2}$: $S_L^{\omega_0}=S_C^{\omega_0}=-\tfrac12$, $S_R^{\omega_0}=0$. $Q=R^{-1}L^{1/2}C^{-1/2}$: $S_R^Q=-1$, $S_L^Q=+\tfrac12$, $S_C^Q=-\tfrac12$.
>
> **Tolerances.**
>
> | Quantity | Worst case | RSS |
> |---|---|---|
> | $\omega_0$ | $\tfrac12(5)+\tfrac12(10)=7.5\,\%$ | $\sqrt{2.5^2+5^2}=5.6\,\%$ |
> | $Q$ | $1+2.5+5=8.5\,\%$ | $\sqrt{1^2+2.5^2+5^2}=5.7\,\%$ |
>
> The square root in $\omega_0$ halves the effect of $L$ and $C$ errors: a nice property, but a tuned circuit with $Q=31.6$ has a $-3\,$dB bandwidth of only $\omega_0/Q\approx3.2\,\%$ of $\omega_0$. A $5.6\,\%$ spread in $\omega_0$ is larger than the passband itself, so such a filter must be trimmed or built with tighter parts.

> [!example] Worked Example 4.49 — The optimum is insensitive: second-order sensitivity
> At an optimum, $\partial y/\partial x=0$, so every first-order sensitivity vanishes. The change is then governed by the second derivative. For the maximum-power problem of Example 4.42, let $r=R_L/R_s$:
> $$
> \frac{P}{P_{\max}}=\frac{4r}{(1+r)^2}.
> $$
> At $r=1$ the first derivative is $0$. The second derivative of $4r(1+r)^{-2}$ is $\dfrac{8(r-2)}{(1+r)^4}$, which equals $-\tfrac12$ at $r=1$, so by Taylor's theorem
> $$
> \frac{P}{P_{\max}}\approx1-\frac14(r-1)^2 .
> $$
> **Numbers.** A $20\,\%$ mismatch ($r=1.2$) is predicted to cost $0.25\times0.04=1\,\%$ of the power; exactly, $4(1.2)/2.2^2=0.9917$, a loss of $0.83\,\%$. ✓
>
> This is a general and useful fact: **near an optimum, errors in the design variable cost only quadratically**. Precise matching is rarely worth paying for when the goal is maximum power; it matters much more when the goal is minimum *reflection*, which depends on mismatch linearly.

> [!example] Worked Example 4.50 — Negative feedback desensitizes gain
> A non-inverting amplifier has closed-loop gain $G=\dfrac{A}{1+A\beta}$, where $A$ is the op-amp's open-loop gain and $\beta$ the feedback fraction set by resistors. With $A=10^5$ and $\beta=0.01$, find $S_A^G$ and $S_\beta^G$.
>
> **Sensitivity to $A$.** $G=A\cdot(1+A\beta)^{-1}$, so by the product, power and sum rules
> $$
> S_A^G=1-S_A^{1+A\beta}=1-\frac{A\beta}{1+A\beta}=\frac{1}{1+A\beta}.
> $$
> **Sensitivity to $\beta$.** $S_\beta^G=-S_\beta^{1+A\beta}=-\dfrac{A\beta}{1+A\beta}$.
>
> **Numbers.** Loop gain $A\beta=1000$: $G=99.90$, $S_A^G=1/1001\approx0.001$, $S_\beta^G\approx-0.999$.
>
> **Interpretation.** A $10\,\%$ change in the op-amp gain changes the closed-loop gain by only about $0.01\,\%$: feedback divides the sensitivity by $1+A\beta$. Even a $50\,\%$ drop ($A=5\times10^4$, a large change for which the linear estimate predicts $0.05\,\%$) gives $G=99.80$, a change of only $0.10\,\%$. The price is that the gain now depends almost entirely on $\beta$, i.e. on two resistors, which are cheap to make accurate. This is the central reason for negative feedback in amplifier design.

> [!example] Worked Example 4.51 — Sensitivity of a strain-gauge bridge
> A Wheatstone bridge excited by $V=5\,$V has one active strain gauge $R=R_0(1+x)$ and three fixed resistors $R_0$. Its output is
> $$
> v_o(x)=V\left(\frac{R}{R+R_0}-\frac12\right)=V\left(\frac{1+x}{2+x}-\frac12\right)=\frac{V}{2}\cdot\frac{x}{2+x}.
> $$
> Find the absolute sensitivity $dv_o/dx$ at balance and the output for a strain of $10^{-3}$ with gauge factor $2$ (so $x=2\times10^{-3}$).
>
> $\dfrac{dv_o}{dx}=\dfrac V2\cdot\dfrac{(2+x)-x}{(2+x)^2}=\dfrac{V}{(2+x)^2}$, which is $V/4=1.25\,$V per unit $x$ at $x=0$. Hence $v_o\approx1.25\times2\times10^{-3}=2.50\,$mV.
>
> **Exact.** $v_o=2.5\times\dfrac{0.002}{2.002}=2.4975\,$mV. The nonlinearity error is $-0.1\,\%$, about $x/2$ in relative terms, which is why bridges are operated near balance (and why half- and full-bridge arrangements, which cancel this term, are preferred).
>
> Here the *relative* sensitivity $S_x^{v_o}$ is not useful, because the nominal output is zero; absolute sensitivity (V per unit strain) is the meaningful figure of merit. Choose the measure that matches the specification.

> [!tip] Engineering Note — temperature coefficients are sensitivities
> A resistor's temperature coefficient of resistance (TCR) is $\alpha=\dfrac1R\dfrac{dR}{dT}$, a relative rate per kelvin, often quoted in ppm/K. Copper has $\alpha\approx3.9\times10^{-3}\,\text{K}^{-1}$; precision metal-film resistors reach $\pm5\,$ppm/K. By the chain rule of sensitivities, the temperature coefficient of any performance measure is $\dfrac1y\dfrac{dy}{dT}=\sum_kS_{x_k}^{y}\,\alpha_k$. For the RC filter of Example 4.47, $f_c$ drifts at $-(\alpha_R+\alpha_C)$; choosing a capacitor whose TC is equal and opposite to the resistor's ($\alpha_C=-\alpha_R$) makes the cutoff frequency first-order temperature-independent. Larger tolerance studies with many parameters and nonlinear effects use Monte Carlo simulation ([Ch. 24](../part-7-numerical-computational/24-computational-mathematics.md)), but the sensitivities computed here tell you which parameters matter before you start.

> [!warning] Common Pitfall — linear sensitivity has a range of validity
> First-order sensitivity assumes small relative changes. With $\pm10\,\%$ or $\pm20\,\%$ tolerances, or near an optimum where the first derivative vanishes, include second-order terms or evaluate the extremes exactly (Examples 4.47 and 4.49). The RSS formula also assumes **independent** errors; components from the same reel, or resistors on one thin-film network, are correlated, and their errors may add (worst case) or cancel (Example 4.46) far more than RSS predicts. Finally, do not confuse absolute and relative sensitivity: $\partial y/\partial x$ has units and depends on the units chosen, $S_x^y$ does not.

## Summary

| Topic | Key result |
|---|---|
| Limit | $\lim_{x\to a}f(x)=L$: $\forall\varepsilon>0\ \exists\delta>0$, $0<\lvert x-a\rvert<\delta\Rightarrow\lvert f(x)-L\rvert<\varepsilon$ |
| Fundamental limits | $\lim_{x\to0}\frac{\sin x}{x}=1$, $\lim_{x\to0}\frac{1-\cos x}{x^2}=\frac12$, $\lim_{h\to0}\frac{e^h-1}{h}=1$, $\lim_{n\to\infty}(1+\frac xn)^n=e^x$ |
| Continuity | $\lim_{x\to a}f(x)=f(a)$; IVT (roots by sign change), EVT (extrema on $[a,b]$) |
| Derivative | $f'(x)=\lim_{h\to0}\frac{f(x+h)-f(x)}{h}$; differentiable $\Rightarrow$ continuous |
| Product, quotient | $(fg)'=f'g+fg'$, $(f/g)'=(f'g-fg')/g^2$ |
| Complex exponential | $\frac{d}{dt}e^{st}=se^{st}$; $\frac{d}{dt}\leftrightarrow j\omega$ for phasors |
| Mean value theorem | $f'(c)=\frac{f(b)-f(a)}{b-a}$ for some $c\in(a,b)$ |
| L'Hôpital | $0/0$ or $\infty/\infty$: $\lim f/g=\lim f'/g'$ |
| Chain rule | $\frac{d}{dx}f(g(x))=f'(g(x))g'(x)$, $\frac{dy}{dx}=\frac{dy}{du}\frac{du}{dx}$ |
| Logarithmic derivative | $\frac{y'}{y}=(\ln y)'$; for $y=\prod f_k^{a_k}$, $\frac{y'}{y}=\sum a_k\frac{f_k'}{f_k}$ |
| Inverse function | $\frac{dx}{dy}=1/\frac{dy}{dx}$ |
| Parametric | $\frac{dy}{dx}=\frac{\dot y}{\dot x}$ |
| Implicit | $F(x,y)=0\Rightarrow y'=-F_x/F_y$ ($F_y\neq0$) |
| Partial derivative | $f_x=\lim_{h\to0}\frac{f(x+h,y)-f(x,y)}{h}$; $f_{xy}=f_{yx}$ if continuous |
| Multivariable chain rule | $\frac{dz}{dt}=f_x\dot x+f_y\dot y$; $\Delta f\approx f_x\Delta x+f_y\Delta y$ |
| Leibniz rule | $(fg)^{(n)}=\sum_{k=0}^n\binom nk f^{(k)}g^{(n-k)}$ |
| Curvature | $\kappa=\lvert f''\rvert/(1+f'^2)^{3/2}$ |
| Taylor series | $f(x)=\sum_{k=0}^n\frac{f^{(k)}(a)}{k!}(x-a)^k+\frac{f^{(n+1)}(c)}{(n+1)!}(x-a)^{n+1}$ |
| Key series | $e^x=\sum\frac{x^n}{n!}$, $\sin x=x-\frac{x^3}{6}+\cdots$, $\cos x=1-\frac{x^2}{2}+\cdots$, $(1+x)^\alpha=1+\alpha x+\cdots$ |
| Finite differences | central: $\frac{f(x+h)-f(x-h)}{2h}=f'+O(h^2)$; second: $\frac{f(x+h)-2f(x)+f(x-h)}{h^2}=f''+O(h^2)$ |
| Linearization | $f(x)\approx f(a)+f'(a)(x-a)$; diode $g_d=I_D/(nV_T)$; MOSFET $g_m=2I_D/V_{ov}$ |
| Equilibrium stability | $\dot x=f(x)$, $f(x^*)=0$: stable if $f'(x^*)<0$, $\tau=-1/f'(x^*)$ |
| Newton's method | $x_{n+1}=x_n-f(x_n)/f'(x_n)$; $e_{n+1}\approx\frac{f''}{2f'}e_n^2$ |
| Rates | $i=\dot q$, $i=\frac{d(Cv)}{dt}$, $v=L\frac{di}{dt}$, $\mathcal E=-\dot\Phi$, max $\lvert\dot v\rvert$ of a sine $=2\pi fV_m$ |
| Optimization | $f'(c)=0$; $f''(c)<0$ max, $>0$ min; check endpoints; max power at $R_L=R_s$, $P_{\max}=V^2/(4R_s)$ |
| Relative sensitivity | $S_x^y=\frac{x}{y}\frac{\partial y}{\partial x}=\frac{\partial\ln y}{\partial\ln x}$; for $y=cx^n$, $S=n$; $S^{uv}=S^u+S^v$ |
| Tolerance | worst case $\sum\lvert S_k\rvert t_k$; RSS $\sqrt{\sum S_k^2\sigma_k^2}$; feedback $S_A^G=1/(1+A\beta)$ |

## Further reading

- J. Stewart, *Calculus: Early Transcendentals*, Cengage — clear, thorough treatment of limits, derivatives, Taylor series and applications, with many geometric explanations.
- E. Kreyszig, *Advanced Engineering Mathematics*, Wiley — the standard engineering reference; useful for the bridge from this chapter to multivariable calculus and ODEs.
- M. Spivak, *Calculus*, Cambridge University Press — rigorous proofs of the theorems sketched here (limits, MVT, Taylor's theorem), for the reader who wants every detail.
- A. S. Sedra and K. C. Smith, *Microelectronic Circuits*, Oxford University Press — small-signal modelling and distortion of diodes and transistors, the main engineering consumer of linearization and Taylor expansion.
