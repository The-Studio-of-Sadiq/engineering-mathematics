# Appendix E — Probability, Statistics and Numerical Methods

[← Table of Contents](../00-front/table-of-contents.md)

> This appendix supplies complete proofs of the central results of Part VI (Chapters 19–21) and of the numerical chapters of Part VII (Chapters 22 and 24) that the chapters state with only a sketch. The probability results run from the elementary (Bayes' theorem, linearity of expectation) through the limit theorems (the weak law of large numbers and the central limit theorem) to the second-order theory of random processes (Wiener–Khinchin, filtering of noise). The statistics results are the unbiasedness of the sample variance, least squares and the Gauss–Markov theorem, and the GUM law of propagation of uncertainty. The numerical results are the convergence theorems for bisection, Newton–Raphson, the trapezoidal and Simpson rules, Euler's method and Monte Carlo integration.
>
> **Prerequisites.** Calculus with Taylor's theorem and its remainder (Ch. 4–6), matrix algebra (Ch. 8), Fourier transforms (Ch. 16), and the probability vocabulary of Ch. 19.
>
> **Level of rigour.** Every hypothesis is stated precisely, and every step is justified from those hypotheses. A few deep tools of measure-theoretic probability and analysis are *quoted* with a reference rather than proved: the linearity and monotonicity of the Lebesgue integral, Fubini's theorem, the dominated convergence theorem, Lévy's continuity theorem, and Fourier inversion. Each is named where it is used. Random variables are real-valued and defined on a probability space $(\Omega,\mathcal F,P)$ (Section 19.1). "i.i.d." means independent and identically distributed.

## E.1 The law of total probability and Bayes' theorem

> [!definition] Conditional probability and partitions
> For events $A,B$ with $P(B)>0$, the conditional probability of $A$ given $B$ is
> $$
> P(A\mid B)=\frac{P(A\cap B)}{P(B)} .
> $$
> A finite or countable family of events $B_1,B_2,\dots$ is a **partition** of $\Omega$ if the $B_i$ are pairwise disjoint ($B_i\cap B_k=\varnothing$ for $i\ne k$) and $\bigcup_iB_i=\Omega$.

**Statement**

> [!theorem] Law of total probability and Bayes' theorem
> Let $B_1,B_2,\dots$ be a finite or countable partition of $\Omega$ with $P(B_i)>0$ for every $i$, and let $A$ be any event.
>
> 1. **Total probability:**
> $$
> P(A)=\sum_i P(A\mid B_i)\,P(B_i).
> $$
> 2. **Bayes' theorem:** if in addition $P(A)>0$, then for every $k$
> $$
> P(B_k\mid A)=\frac{P(A\mid B_k)\,P(B_k)}{\sum_i P(A\mid B_i)\,P(B_i)} .
> $$

**Proof**

*Step 1 (decompose $A$).* Because the $B_i$ cover $\Omega$,

$$
A=A\cap\Omega=A\cap\Big(\bigcup_iB_i\Big)=\bigcup_i\,(A\cap B_i),
$$

by the distributive law of intersection over union.

*Step 2 (the pieces are disjoint).* For $i\ne k$, $(A\cap B_i)\cap(A\cap B_k)\subseteq B_i\cap B_k=\varnothing$.

*Step 3 (additivity).* By the countable additivity axiom (Section 19.1), applied to the disjoint family of Step 2,

$$
P(A)=\sum_iP(A\cap B_i).
$$

*Step 4 (multiplication rule).* Since $P(B_i)>0$, the definition of conditional probability rearranges to $P(A\cap B_i)=P(A\mid B_i)P(B_i)$. Substituting into Step 3 proves part 1.

*Step 5 (Bayes).* Since $P(A)>0$ and $P(B_k)>0$, the definition gives two expressions for the same joint probability:

$$
P(B_k\mid A)\,P(A)=P(A\cap B_k)=P(A\mid B_k)\,P(B_k).
$$

Divide by $P(A)$ and replace $P(A)$ by the sum from part 1. ∎

**Remarks**

- *Why the hypotheses matter.* If the $B_i$ overlap, Step 2 fails and the sum in part 1 counts the overlap twice. For example, with $B_1=B_2=\Omega$ the formula would give $P(A)=2P(A)$. If the $B_i$ do not cover $\Omega$, the sum misses $A\setminus\bigcup_iB_i$. The condition $P(B_i)>0$ is needed only so that $P(A\mid B_i)$ is defined. Cells with $P(B_i)=0$ may simply be dropped, because they contribute $P(A\cap B_i)=0$.
- *Base-rate effect.* Suppose a fault detector flags a faulty unit with probability $0.99$, flags a healthy unit with probability $0.05$, and $1\%$ of units are faulty. Bayes' theorem gives

$$
P(\text{faulty}\mid\text{flag})=\frac{0.99\cdot0.01}{0.99\cdot0.01+0.05\cdot0.99}=\frac16\approx0.167 .
$$

So five out of six alarms are false, because the small prior is weighted against a comparatively large false-alarm rate. The chapter treats this in [Section 19.3](../part-6-probability-statistics/19-probability.md).
- *Continuous form.* If $X$ has density $f_X$ and the conditional densities exist, the same argument with integrals in place of sums gives $f_{X\mid A}(x)=P(A\mid X=x)f_X(x)/\int P(A\mid X=u)f_X(u)\,du$. Strictly, conditioning on the null event $\{X=x\}$ needs the measure-theoretic definition of conditional probability. The formula used in [Section 19.7](../part-6-probability-statistics/19-probability.md) is the limit of the discrete one over shrinking cells $\{x<X\le x+\Delta\}$.
- *Where used.* Total probability appears in [Section 19.2](../part-6-probability-statistics/19-probability.md). Bayes' theorem underlies the binary channel and MAP/ML detection in [Section 19.3](../part-6-probability-statistics/19-probability.md), and likelihood ratio tests in [Section 20.8](../part-6-probability-statistics/20-statistical-analysis.md).

## E.2 Linearity of expectation and the variance of a sum

> [!definition] Expectation, variance, covariance
> A random variable $X$ is **integrable** if $E\lvert X\rvert<\infty$. Its expectation is $E[X]=\int_\Omega X\,dP$. This equals $\sum_xx\,p_X(x)$ for a pmf and $\int xf_X(x)\,dx$ for a pdf (Section 19.9). If $E[X^2]<\infty$, the **variance** is $\operatorname{Var}(X)=E[(X-\mu_X)^2]$ with $\mu_X=E[X]$. For $X,Y$ with finite second moments, the **covariance** is $\operatorname{Cov}(X,Y)=E[(X-\mu_X)(Y-\mu_Y)]$.

**Statement**

> [!theorem] Linearity of expectation; variance of a sum
> 1. If $X$ and $Y$ are integrable and $a,b\in\mathbb R$, then $aX+bY$ is integrable and
> $$
> E[aX+bY]=a\,E[X]+b\,E[Y].
> $$
> No independence is assumed.
> 2. If $E[X^2]<\infty$ and $E[Y^2]<\infty$, then $\operatorname{Cov}(X,Y)$ is finite and
> $$
> \operatorname{Var}(X+Y)=\operatorname{Var}(X)+\operatorname{Var}(Y)+2\operatorname{Cov}(X,Y).
> $$
> 3. More generally, for a random vector $\mathbf X=(X_1,\dots,X_n)^\mathsf T$ with finite second moments, covariance matrix $\boldsymbol\Sigma$ with entries $\Sigma_{ik}=\operatorname{Cov}(X_i,X_k)$, and $\mathbf a\in\mathbb R^n$,
> $$
> \operatorname{Var}(\mathbf a^\mathsf T\mathbf X)=\sum_{i=1}^n\sum_{k=1}^na_ia_k\Sigma_{ik}=\mathbf a^\mathsf T\boldsymbol\Sigma\,\mathbf a ,
> $$
> and $\boldsymbol\Sigma$ is symmetric positive semidefinite.

**Proof**

*Part 1.* Integrability follows from the triangle inequality, $\lvert aX+bY\rvert\le\lvert a\rvert\lvert X\rvert+\lvert b\rvert\lvert Y\rvert$, together with monotonicity of the integral. The identity is the linearity of the Lebesgue integral $\int_\Omega(\cdot)\,dP$, which we quote. The two cases used in this book can also be checked directly.

- *Joint pmf $p(x,y)$.* By the law of the unconscious statistician (Section 19.9),

$$
E[aX+bY]=\sum_{x,y}(ax+by)\,p(x,y)=a\sum_xx\sum_yp(x,y)+b\sum_yy\sum_xp(x,y)=a\sum_xx\,p_X(x)+b\sum_yy\,p_Y(y).
$$

The rearrangement is legitimate because the double series converges absolutely: $\sum\lvert ax+by\rvert\,p\le\lvert a\rvert E\lvert X\rvert+\lvert b\rvert E\lvert Y\rvert<\infty$.
- *Joint pdf $f(x,y)$.* The same computation with $\iint$ in place of $\sum\sum$ works, using Fubini's theorem (justified by the same absolute bound) and the marginalisation $\int f(x,y)\,dy=f_X(x)$.

Only the marginals $p_X$ and $p_Y$ survive the computation. The dependence between $X$ and $Y$, which is encoded in $p(x,y)$, has been summed out. That is why independence is not needed.

*Part 2.* First, finiteness. From $(\lvert x\rvert-\lvert y\rvert)^2\ge0$ we get $\lvert xy\rvert\le\frac12(x^2+y^2)$, so $E\lvert XY\rvert\le\frac12(E[X^2]+E[Y^2])<\infty$. Likewise $E\lvert X\rvert\le\frac12(1+E[X^2])<\infty$. All the expectations below therefore exist. Put $\tilde X=X-\mu_X$ and $\tilde Y=Y-\mu_Y$. By part 1, $E[X+Y]=\mu_X+\mu_Y$, so $X+Y-E[X+Y]=\tilde X+\tilde Y$ and

$$
\operatorname{Var}(X+Y)=E\big[(\tilde X+\tilde Y)^2\big]=E\big[\tilde X^2+\tilde Y^2+2\tilde X\tilde Y\big]=E[\tilde X^2]+E[\tilde Y^2]+2E[\tilde X\tilde Y],
$$

where the last equality is part 1 again. The three terms are $\operatorname{Var}X$, $\operatorname{Var}Y$ and $\operatorname{Cov}(X,Y)$.

*Part 3.* With $\tilde X_i=X_i-\mu_i$, we have $\mathbf a^\mathsf T\mathbf X-E[\mathbf a^\mathsf T\mathbf X]=\sum_ia_i\tilde X_i$. Squaring the finite sum and applying part 1 term by term gives

$$
E\Big[\Big(\sum_ia_i\tilde X_i\Big)^2\Big]=\sum_i\sum_ka_ia_kE[\tilde X_i\tilde X_k]=\mathbf a^\mathsf T\boldsymbol\Sigma\,\mathbf a .
$$

Symmetry of $\boldsymbol\Sigma$ is clear from the definition. Positive semidefiniteness holds because $\mathbf a^\mathsf T\boldsymbol\Sigma\mathbf a$ is a variance, the expectation of a square, and so is $\ge0$ for every $\mathbf a$. ∎

**Remarks**

- *Independent or uncorrelated terms.* If the $X_i$ are pairwise **uncorrelated** ($\Sigma_{ik}=0$ for $i\ne k$), part 3 reduces to $\operatorname{Var}(\sum a_iX_i)=\sum a_i^2\operatorname{Var}X_i$. Independence implies uncorrelatedness (Section 19.9) but is not needed. This is the "quadrature addition" of tolerances and noise powers.
- *Counterexample to "variances always add".* Take $Y=-X$. Then $\operatorname{Var}(X+Y)=0$, while $\operatorname{Var}X+\operatorname{Var}Y=2\operatorname{Var}X$. The covariance term $-2\operatorname{Var}X$ accounts for the difference exactly. Differential signalling exploits this. Common-mode noise is positively correlated on the two wires, so it cancels in the difference.
- *Finite moments are essential.* If $X$ is Cauchy, $E\lvert X\rvert=\infty$, and $E[X]+E[-X]$ is the meaningless $\infty-\infty$, even though $X+(-X)=0$.
- *Where used.* [Section 19.9](../part-6-probability-statistics/19-probability.md) (sums, tolerance stack-up), [Section 20.3](../part-6-probability-statistics/20-statistical-analysis.md) (covariance matrix), and every later section of this appendix.

## E.3 Markov and Chebyshev inequalities

**Statement**

> [!theorem] Markov and Chebyshev inequalities
> 1. **Markov.** If $X\ge0$ and $a>0$, then
> $$
> P(X\ge a)\le\frac{E[X]}{a}.
> $$
> 2. **Chebyshev.** If $E[X^2]<\infty$, $\mu=E[X]$, $\sigma^2=\operatorname{Var}X$ and $\varepsilon>0$, then
> $$
> P(\lvert X-\mu\rvert\ge\varepsilon)\le\frac{\sigma^2}{\varepsilon^2},\qquad\text{equivalently}\qquad P(\lvert X-\mu\rvert\ge k\sigma)\le\frac1{k^2}\quad(k>0,\ \sigma>0).
> $$

**Proof**

*Markov.* Let $I=\mathbf 1_{\{X\ge a\}}$ be the indicator of the event $\{X\ge a\}$. We claim that $aI\le X$ at every outcome $\omega$. If $X(\omega)\ge a$, then $aI(\omega)=a\le X(\omega)$. Otherwise $aI(\omega)=0\le X(\omega)$, because $X\ge0$. Expectation is monotone: if $U\le V$, then $E[V]-E[U]=E[V-U]\ge0$ by E.2 and the non-negativity of the integral of a non-negative function. Hence $E[aI]\le E[X]$. But $E[aI]=a\,P(X\ge a)$. Divide by $a>0$. If $E[X]=\infty$, there is nothing to prove.

*Chebyshev.* The variable $Y=(X-\mu)^2$ is non-negative with $E[Y]=\sigma^2$. The events $\{\lvert X-\mu\rvert\ge\varepsilon\}$ and $\{Y\ge\varepsilon^2\}$ are identical, because $t\mapsto t^2$ is increasing on $t\ge0$. Markov with $a=\varepsilon^2$ gives $P(Y\ge\varepsilon^2)\le\sigma^2/\varepsilon^2$. Setting $\varepsilon=k\sigma$ gives the second form. ∎

```desmos-graph
left=-0.5; right=6; bottom=-0.5; top=5
---
y=x|blue|0<x<6
y=0|red|0<x<2
y=2|red|2<x<6
(2,2)|red|label:a=2
```

*The pointwise inequality behind Markov's bound for $a=2$: the step $a\,\mathbf 1_{\{x\ge a\}}$ (red) never exceeds $x$ (blue) on $x\ge0$. Taking expectations of both sides gives $aP(X\ge a)\le E[X]$.*

**Remarks**

- *Sharpness.* Both inequalities are attained, so they cannot be improved without more information about the distribution. For Markov, let $X=a$ with probability $p$ and $X=0$ otherwise. Then $E[X]/a=p=P(X\ge a)$. For Chebyshev with $k\ge1$, take the three-point law $P(X=\pm k\sigma)=\frac1{2k^2}$, $P(X=0)=1-\frac1{k^2}$. It has mean $0$, variance $\sigma^2$ and $P(\lvert X\rvert\ge k\sigma)=1/k^2$.
- *Non-negativity matters for Markov.* For $X=\pm1$ with equal probability, $E[X]=0$ but $P(X\ge1)=\frac12>0$. The proof breaks at $aI\le X$ when $X(\omega)=-1$.
- *Generalisations.* Applying Markov to $\lvert X-\mu\rvert^r$ gives $P(\lvert X-\mu\rvert\ge\varepsilon)\le E\lvert X-\mu\rvert^r/\varepsilon^r$. Applying it to $e^{sX}$ gives the Chernoff bound $P(X\ge a)\le e^{-sa}M_X(s)$ for $s>0$, the basis of the exponential error bounds of coding theory.
- *Where used.* Distribution-free yield guarantees ([Section 19.9](../part-6-probability-statistics/19-probability.md), Example 19.43), the weak law in E.4, and the non-asymptotic Monte Carlo error bound in E.14.

## E.4 The weak law of large numbers

> [!definition] Convergence in probability and in distribution
> $Z_n\to c$ **in probability** if $P(\lvert Z_n-c\rvert\ge\varepsilon)\to0$ as $n\to\infty$ for every $\varepsilon>0$. $Z_n\to Z$ **in distribution** if $F_{Z_n}(z)\to F_Z(z)$ at every point $z$ where $F_Z$ is continuous.

**Statement**

> [!theorem] Weak law of large numbers
> Let $\bar X_n=\frac1n\sum_{i=1}^nX_i$.
>
> 1. **($L^2$ version.)** If $X_1,X_2,\dots$ are pairwise uncorrelated with common mean $\mu$ and variances $\operatorname{Var}X_i\le\sigma^2<\infty$, then for every $\varepsilon>0$
> $$
> P\big(\lvert\bar X_n-\mu\rvert\ge\varepsilon\big)\le\frac{\sigma^2}{n\varepsilon^2}\ \longrightarrow\ 0 .
> $$
> 2. **(Khinchin's version.)** If $X_1,X_2,\dots$ are i.i.d. with $E\lvert X_1\rvert<\infty$ and $\mu=E[X_1]$, then $\bar X_n\to\mu$ in probability. No variance is assumed.

**Proof**

*Part 1.* By linearity (E.2), $E[\bar X_n]=\frac1n\sum E[X_i]=\mu$. Because the $X_i$ are pairwise uncorrelated, E.2 part 3 with $a_i=1/n$ gives

$$
\operatorname{Var}(\bar X_n)=\frac1{n^2}\sum_{i=1}^n\operatorname{Var}X_i\le\frac{n\sigma^2}{n^2}=\frac{\sigma^2}{n}.
$$

Chebyshev's inequality (E.3) applied to $\bar X_n$ gives $P(\lvert\bar X_n-\mu\rvert\ge\varepsilon)\le\operatorname{Var}(\bar X_n)/\varepsilon^2\le\sigma^2/(n\varepsilon^2)$.

*Part 2.* We use characteristic functions $\varphi_X(\omega)=E[e^{j\omega X}]$ (Section 19.9), Lévy's continuity theorem (stated in E.5), the product Lemma E.5a (proved in E.5), and one elementary fact:

- (ii) If $Z_n\to c$ in distribution and $c$ is a constant, then $Z_n\to c$ in probability.

For (ii), note that the limit CDF $\mathbf 1_{\{z\ge c\}}$ is continuous at $c-\varepsilon$ and at $c+\varepsilon/2$. So $P(\lvert Z_n-c\rvert\ge\varepsilon)\le F_{Z_n}(c-\varepsilon)+1-F_{Z_n}(c+\varepsilon/2)\to0+1-1=0$.

*Step 1 (first-order expansion).* Let $\varphi=\varphi_{X_1}$. For real $x$, $\lvert e^{jx}-1\rvert\le\min(\lvert x\rvert,2)$, and hence $\lvert e^{jx}-1-jx\rvert\le\min(x^2/2,\,2\lvert x\rvert)$. The bounds $\lvert x\rvert$ and $x^2/2$ come from writing $e^{jx}-1=\int_0^xje^{ju}du$ and $e^{jx}-1-jx=\int_0^xj(e^{ju}-1)du$ and using $\lvert e^{ju}\rvert=1$. The bounds $2$ and $2\lvert x\rvert$ come from the triangle inequality. Hence

$$
\big\lvert\varphi(\omega)-1-j\omega\mu\big\rvert=\big\lvert E\big[e^{j\omega X_1}-1-j\omega X_1\big]\big\rvert\le\lvert\omega\rvert\,E\Big[\min\Big(\tfrac{\lvert\omega\rvert X_1^2}{2},\,2\lvert X_1\rvert\Big)\Big].
$$

The random variable inside the last expectation tends to $0$ pointwise as $\omega\to0$ and is dominated by $2\lvert X_1\rvert$, which is integrable. By dominated convergence, the expectation tends to $0$. Therefore $\varphi(\omega)=1+j\omega\mu+\omega\,r(\omega)$ with $r(\omega)\to0$ as $\omega\to0$.

*Step 2 (product).* By independence, the characteristic function of a sum is the product of the characteristic functions, so $\varphi_{\bar X_n}(\omega)=\varphi(\omega/n)^n$. Fix $\omega$. Lemma E.5a needs factors of modulus at most $1$. Characteristic functions satisfy $\lvert\varphi\rvert\le1$, and so does $e^{j\omega\mu/n}$. Using Step 1 and $\lvert e^{jx}-1-jx\rvert\le x^2/2$,

$$
\big\lvert\varphi(\omega/n)-e^{j\omega\mu/n}\big\rvert\le\frac{\lvert\omega\rvert}{n}\,\lvert r(\omega/n)\rvert+\frac{\omega^2\mu^2}{2n^2}.
$$

Lemma E.5a with $n$ equal factors then gives

$$
\big\lvert\varphi(\omega/n)^n-e^{j\omega\mu}\big\rvert\le n\,\big\lvert\varphi(\omega/n)-e^{j\omega\mu/n}\big\rvert\le\lvert\omega\rvert\,\lvert r(\omega/n)\rvert+\frac{\omega^2\mu^2}{2n}\ \longrightarrow\ 0 .
$$

*Step 3 (conclusion).* $e^{j\omega\mu}$ is the characteristic function of the constant $\mu$. By Lévy's theorem $\bar X_n\to\mu$ in distribution, and by (ii) in probability. ∎

**Remarks**

- *What the law says and does not say.* Part 1 gives a *rate*, $P(\lvert\bar X_n-\mu\rvert\ge\varepsilon)=O(1/n)$, with an explicit constant. It justifies the relative-frequency interpretation of [Section 19.1](../part-6-probability-statistics/19-probability.md). With $X_i=\mathbf 1_A$ on independent trials, $\bar X_n=N_A/n\to P(A)$. It says nothing about a *single* infinite run of experiments converging. That is the strong law ($\bar X_n\to\mu$ with probability one), whose proof is harder and which is quoted in [Section 24.6](../part-7-numerical-computational/24-computational-mathematics.md).
- *Counterexample without a mean.* If the $X_i$ are i.i.d. standard Cauchy, then $\varphi(\omega)=e^{-\lvert\omega\rvert}$ and $\varphi(\omega/n)^n=e^{-\lvert\omega\rvert}$. So $\bar X_n$ is again standard Cauchy for every $n$ and does not concentrate at all. Averaging does not help against impulsive noise of this kind. Median filters do.
- *Correlation matters.* If $X_i=X_1$ for all $i$ (perfect correlation), then $\bar X_n=X_1$ and nothing converges. Part 1 needs the variance of the average to vanish. The same principle, with correlation that decays in time, gives the mean-ergodic theorem of [Section 21.7](../part-6-probability-statistics/21-random-processes.md).
- *Where used.* Consistency of estimators ([Section 20.6](../part-6-probability-statistics/20-statistical-analysis.md)), Monte Carlo (E.14), ergodicity ([Section 21.7](../part-6-probability-statistics/21-random-processes.md)).

## E.5 The central limit theorem

**Statement**

> [!theorem] Central limit theorem (Lindeberg–Lévy)
> Let $X_1,X_2,\dots$ be i.i.d. with mean $\mu$ and finite variance $\sigma^2>0$. Put
> $$
> Z_n=\frac{\sum_{i=1}^nX_i-n\mu}{\sigma\sqrt n}=\frac{\bar X_n-\mu}{\sigma/\sqrt n}.
> $$
> Then $Z_n\to\mathcal N(0,1)$ in distribution, i.e. $P(Z_n\le z)\to\Phi(z)$ for every $z\in\mathbb R$.

The proof rests on one quoted theorem and two lemmas.

> [!theorem] Lévy's continuity theorem (quoted)
> Let $Z_n$ have characteristic functions $\varphi_n$. If $\varphi_n(\omega)\to\psi(\omega)$ for every $\omega\in\mathbb R$, and $\psi$ is continuous at $\omega=0$, then $\psi$ is the characteristic function of some random variable $Z$, and $Z_n\to Z$ in distribution.

For a proof, see Billingsley (*Probability and Measure*, §26) or Durrett (*Probability: Theory and Examples*, §3.3). It is the probabilistic counterpart of the uniqueness of the Fourier transform (Ch. 16). The characteristic function is the Fourier transform of the pdf, with the sign of the exponent reversed.

> [!theorem] Lemma E.5a (products of numbers in the unit disc)
> If $a_1,\dots,a_n,b_1,\dots,b_n\in\mathbb C$ all have modulus at most $1$, then
> $$
> \Big\lvert\prod_{i=1}^na_i-\prod_{i=1}^nb_i\Big\rvert\le\sum_{i=1}^n\lvert a_i-b_i\rvert .
> $$

*Proof.* Use induction on $n$. The case $n=1$ is trivial. For the inductive step, write $A_m=\prod_{i\le m}a_i$ and $B_m=\prod_{i\le m}b_i$, both of modulus $\le1$. Then

$$
A_n-B_n=a_n(A_{n-1}-B_{n-1})+(a_n-b_n)B_{n-1},
$$

so $\lvert A_n-B_n\rvert\le\lvert A_{n-1}-B_{n-1}\rvert+\lvert a_n-b_n\rvert$, which together with the induction hypothesis proves the claim. ∎

> [!theorem] Lemma E.5b (second-order expansion of a characteristic function)
> If $E[Y]=0$ and $E[Y^2]=1$, then $\varphi_Y(\omega)=1-\frac{\omega^2}{2}+\omega^2r(\omega)$ with $r(\omega)\to0$ as $\omega\to0$.

*Proof.* For real $x$,

$$
\Big\lvert e^{jx}-1-jx+\frac{x^2}{2}\Big\rvert\le\min\Big(\frac{\lvert x\rvert^3}{6},\,x^2\Big).
$$

For the first bound, note that $e^{jx}-1-jx+\frac{x^2}2=\int_0^xj\big(e^{ju}-1-ju\big)\,du$. The integrand has modulus at most $u^2/2$ (E.4, Step 1), which integrates to $\lvert x\rvert^3/6$. For the second bound, $\lvert e^{jx}-1-jx\rvert\le x^2/2$ plus $\lvert x^2/2\rvert$ gives $x^2$. Now set $x=\omega Y$ and take expectations. Using $E[Y]=0$ and $E[Y^2]=1$,

$$
\Big\lvert\varphi_Y(\omega)-1+\frac{\omega^2}{2}\Big\rvert\le E\Big[\min\Big(\frac{\lvert\omega\rvert^3\lvert Y\rvert^3}{6},\,\omega^2Y^2\Big)\Big]=\omega^2\,E\Big[\min\Big(\frac{\lvert\omega\rvert\,\lvert Y\rvert^3}{6},\,Y^2\Big)\Big].
$$

Define $r(\omega)=\big(\varphi_Y(\omega)-1+\omega^2/2\big)/\omega^2$ for $\omega\ne0$. The display says that $\lvert r(\omega)\rvert$ is at most the last expectation. The variable inside that expectation tends to $0$ pointwise as $\omega\to0$ and is dominated by the integrable $Y^2$, so by dominated convergence $r(\omega)\to0$. ∎

**Proof of the theorem**

*Step 1 (standardise).* Let $Y_i=(X_i-\mu)/\sigma$. The $Y_i$ are i.i.d. with $E[Y_i]=0$ and $E[Y_i^2]=1$, and $Z_n=n^{-1/2}\sum_{i=1}^nY_i$. Let $\varphi$ be the common characteristic function of the $Y_i$.

*Step 2 (characteristic function of $Z_n$).* By independence, the expectation of a product of functions of independent variables factorises (Section 19.9). Hence

$$
\varphi_{Z_n}(\omega)=E\Big[\prod_{i=1}^ne^{j\omega Y_i/\sqrt n}\Big]=\prod_{i=1}^nE\big[e^{j\omega Y_i/\sqrt n}\big]=\varphi\Big(\frac{\omega}{\sqrt n}\Big)^n .
$$

*Step 3 (compare with a simple product).* Fix $\omega$ and let $n\ge\omega^2/2$, so that $b_n=1-\frac{\omega^2}{2n}\in[0,1]$. Every characteristic function satisfies $\lvert\varphi\rvert\le E\lvert e^{j\omega Y}\rvert=1$. Lemma E.5a with $n$ equal factors and Lemma E.5b give

$$
\Big\lvert\varphi\Big(\frac{\omega}{\sqrt n}\Big)^n-b_n^n\Big\rvert\le n\,\Big\lvert\varphi\Big(\frac{\omega}{\sqrt n}\Big)-1+\frac{\omega^2}{2n}\Big\rvert=n\cdot\frac{\omega^2}{n}\Big\lvert r\Big(\frac{\omega}{\sqrt n}\Big)\Big\rvert=\omega^2\Big\lvert r\Big(\frac{\omega}{\sqrt n}\Big)\Big\rvert\ \longrightarrow\ 0 .
$$

*Step 4 (the limit).* $b_n^n=\big(1-\frac{\omega^2}{2n}\big)^n\to e^{-\omega^2/2}$, the standard limit $(1+c/n)^n\to e^c$ (Ch. 4). Together with Step 3, $\varphi_{Z_n}(\omega)\to e^{-\omega^2/2}$ for every $\omega$.

*Step 5 (identify and conclude).* $e^{-\omega^2/2}$ is continuous and is the characteristic function of $\mathcal N(0,1)$ (Section 19.9). By Lévy's continuity theorem, $Z_n\to\mathcal N(0,1)$ in distribution. Since $\Phi$ is continuous everywhere, $P(Z_n\le z)\to\Phi(z)$ at every $z$. ∎

```desmos-graph
left=-6; right=6; bottom=-0.3; top=1.05
---
y=\frac{\sin(1.7321x)}{1.7321x}|blue
y=\left(\frac{\sin(1.2247x)}{1.2247x}\right)^{2}|green
y=\left(\frac{\sin(0.7746x)}{0.7746x}\right)^{5}|orange
y=e^{-\frac{x^{2}}{2}}|black|dashed
```

*The proof in pictures. The characteristic function $\varphi(\omega)=\sin(\sqrt3\,\omega)/(\sqrt3\,\omega)$ of a unit-variance uniform variable (blue) and $\varphi(\omega/\sqrt n)^n$ for $n=2$ (green) and $n=5$ (orange) close in on $e^{-\omega^2/2}$ (dashed). At $\omega=1.5$ the values are $0.199,\ 0.276,\ 0.307$, and $0.320$ for $n=20$, against $0.325$.*

**Remarks**

- *Finite variance is essential.* For i.i.d. standard Cauchy variables, $\varphi(\omega)=e^{-\lvert\omega\rvert}$ and $\varphi(\omega/\sqrt n)^n=e^{-\sqrt n\lvert\omega\rvert}$. This tends to $0$ for $\omega\ne0$ and to $1$ at $\omega=0$, a limit that is not continuous at $0$, so Lévy's theorem does not apply. Indeed $Z_n$ (with $\mu=0$, $\sigma=1$ inserted formally) is Cauchy with scale $\sqrt n$, and its mass escapes to infinity. The correct normalisation for Cauchy sums is $1/n$, under which $\bar X_n$ stays Cauchy (E.4). Lemma E.5b is exactly where $E[Y^2]<\infty$ is used.
- *Independence is essential.* If $X_i=X_1$ for all $i$, then $Z_n=\sqrt n\,(X_1-\mu)/\sigma$, which never becomes Gaussian unless $X_1$ already is. Weak dependence is allowed in versions for mixing sequences, which matter for correlated noise samples (Section 21.7).
- *Identical distribution can be relaxed.* The Lindeberg condition allows independent, non-identical summands provided no single term dominates the variance. This is why the sum of many small, independent physical contributions (electron motions in a resistor, Section 21.8) is Gaussian.
- *Rate.* The proof gives no rate. The Berry–Esseen theorem, $\sup_z\lvert P(Z_n\le z)-\Phi(z)\rvert\le C\,E\lvert Y\rvert^3/\sqrt n$, needs a third moment ([Section 20.5](../part-6-probability-statistics/20-statistical-analysis.md)).
- *Where used.* Sampling distributions and confidence intervals ([Sections 20.5 and 20.7](../part-6-probability-statistics/20-statistical-analysis.md)), Gaussian noise models ([Section 21.5](../part-6-probability-statistics/21-random-processes.md)), and Monte Carlo error bars (E.14).

## E.6 The sample variance with divisor $n-1$ is unbiased

**Statement**

> [!theorem] Unbiasedness of $S^2$
> Let $X_1,\dots,X_n$ ($n\ge2$) be pairwise uncorrelated with common mean $\mu$ and common finite variance $\sigma^2$. Let $\bar X=\frac1n\sum X_i$ and
> $$
> S^2=\frac1{n-1}\sum_{i=1}^n(X_i-\bar X)^2 .
> $$
> Then $E[S^2]=\sigma^2$.

**Proof**

*Step 1 (algebraic identity).* Write $X_i-\bar X=(X_i-\mu)-(\bar X-\mu)$ and expand the square:

$$
\sum_{i=1}^n(X_i-\bar X)^2=\sum_{i=1}^n(X_i-\mu)^2-2(\bar X-\mu)\sum_{i=1}^n(X_i-\mu)+n(\bar X-\mu)^2 .
$$

Since $\sum_i(X_i-\mu)=n(\bar X-\mu)$, the middle term is $-2n(\bar X-\mu)^2$. Therefore

$$
\sum_{i=1}^n(X_i-\bar X)^2=\sum_{i=1}^n(X_i-\mu)^2-n(\bar X-\mu)^2 .
$$

*Step 2 (expectations).* $E[(X_i-\mu)^2]=\sigma^2$ by definition. $E[\bar X]=\mu$ by linearity, so $E[(\bar X-\mu)^2]=\operatorname{Var}\bar X=\sigma^2/n$ by E.2, part 3, using pairwise uncorrelatedness.

*Step 3 (combine).* By linearity,

$$
E\Big[\sum_{i=1}^n(X_i-\bar X)^2\Big]=n\sigma^2-n\cdot\frac{\sigma^2}{n}=(n-1)\sigma^2,
$$

and dividing by $n-1$ gives $E[S^2]=\sigma^2$. ∎

**Remarks**

- *The divisor $n$ is biased.* $\hat\sigma^2=\frac1n\sum(X_i-\bar X)^2$ has $E[\hat\sigma^2]=\frac{n-1}{n}\sigma^2$, an underestimate by $\sigma^2/n$. A simulation with $n=5$ standard Gaussian samples ($2\times10^5$ replications) gives mean values $0.999$ for $S^2$ and $0.799$ for $\hat\sigma^2$, against the predicted $1$ and $0.8$. The bias comes from Step 1. The data sit closer to their own mean $\bar X$ than to $\mu$, and the deficit is exactly $n(\bar X-\mu)^2$, with expectation $\sigma^2$.
- *Counterexample: correlated samples.* If the $X_i$ are equicorrelated, $\operatorname{Cov}(X_i,X_k)=\rho\sigma^2$ for $i\ne k$, then $\operatorname{Var}\bar X=\sigma^2\big(1+(n-1)\rho\big)/n$. Step 3 becomes $n\sigma^2-\sigma^2\big(1+(n-1)\rho\big)=(n-1)(1-\rho)\sigma^2$, so $E[S^2]=(1-\rho)\sigma^2$. Positively correlated samples, such as noise sampled faster than its correlation time ([Section 21.3](../part-6-probability-statistics/21-random-processes.md)), make $S^2$ underestimate the variance. Slow drift is invisible to it.
- *$S$ is biased for $\sigma$.* $\sqrt{\cdot}$ is strictly concave, so Jensen's inequality gives $E[S]<\sqrt{E[S^2]}=\sigma$ unless $S$ is almost surely constant.
- *Where used.* [Section 20.2](../part-6-probability-statistics/20-statistical-analysis.md) ("Why divide by $n-1$?"), the $\chi^2_{n-1}$ sampling distribution in [Section 20.5](../part-6-probability-statistics/20-statistical-analysis.md), Type A uncertainty in [Section 20.10](../part-6-probability-statistics/20-statistical-analysis.md), and the Monte Carlo standard error in [Section 24.6](../part-7-numerical-computational/24-computational-mathematics.md). The regression analogue, divisor $n-p$, is proved in E.7.

## E.7 Least-squares regression and the Gauss–Markov theorem

> [!definition] Linear model
> The **linear model** is $\mathbf y=\mathbf X\boldsymbol\beta+\boldsymbol\varepsilon$. Here $\mathbf y\in\mathbb R^n$ is observed, the **design matrix** $\mathbf X\in\mathbb R^{n\times p}$ is known and deterministic, $\boldsymbol\beta\in\mathbb R^p$ is unknown, and the errors satisfy $E[\boldsymbol\varepsilon]=\mathbf 0$ and $\operatorname{Cov}(\boldsymbol\varepsilon)=\sigma^2\mathbf I$, i.e. they are uncorrelated with equal variance (*homoscedastic*). Simple linear regression $Y_i=\beta_0+\beta_1x_i+\varepsilon_i$ is the case $p=2$, where row $i$ of $\mathbf X$ is $(1,\ x_i)$.

**Statement**

> [!theorem] Least squares and its sampling properties
> Suppose $\operatorname{rank}\mathbf X=p$ (so $n\ge p$).
>
> 1. $\mathbf X^\mathsf T\mathbf X$ is invertible, and $S(\boldsymbol\beta)=\lVert\mathbf y-\mathbf X\boldsymbol\beta\rVert^2$ has the **unique** minimiser
> $$
> \hat{\boldsymbol\beta}=(\mathbf X^\mathsf T\mathbf X)^{-1}\mathbf X^\mathsf T\mathbf y .
> $$
> The residual $\mathbf e=\mathbf y-\mathbf X\hat{\boldsymbol\beta}$ satisfies $\mathbf X^\mathsf T\mathbf e=\mathbf 0$. For simple regression, the rank condition says that the $x_i$ are not all equal ($S_{xx}>0$), and
> $$
> \hat\beta_1=\frac{S_{xy}}{S_{xx}},\qquad\hat\beta_0=\bar y-\hat\beta_1\bar x,\qquad S_{xx}=\sum(x_i-\bar x)^2,\ \ S_{xy}=\sum(x_i-\bar x)(y_i-\bar y).
> $$
> 2. Under the linear model, $E[\hat{\boldsymbol\beta}]=\boldsymbol\beta$, $\operatorname{Cov}(\hat{\boldsymbol\beta})=\sigma^2(\mathbf X^\mathsf T\mathbf X)^{-1}$, and $E[\mathrm{SSE}]=(n-p)\sigma^2$, where $\mathrm{SSE}=\lVert\mathbf e\rVert^2$. In simple regression, $\operatorname{Var}\hat\beta_1=\sigma^2/S_{xx}$ and $\operatorname{Var}\hat\beta_0=\sigma^2\sum x_i^2/(nS_{xx})$.
> 3. **Gauss–Markov.** Let $\tilde{\boldsymbol\beta}=\mathbf C\mathbf y$, with $\mathbf C\in\mathbb R^{p\times n}$ deterministic, be any linear estimator that is unbiased for every $\boldsymbol\beta$. Then $\operatorname{Cov}(\tilde{\boldsymbol\beta})-\operatorname{Cov}(\hat{\boldsymbol\beta})$ is positive semidefinite. In particular, $\operatorname{Var}(\mathbf a^\mathsf T\tilde{\boldsymbol\beta})\ge\operatorname{Var}(\mathbf a^\mathsf T\hat{\boldsymbol\beta})$ for every $\mathbf a\in\mathbb R^p$. Least squares is the **best linear unbiased estimator** (BLUE).

**Proof**

*Part 1 (invertibility).* If $\mathbf X^\mathsf T\mathbf X\mathbf v=\mathbf 0$, then $0=\mathbf v^\mathsf T\mathbf X^\mathsf T\mathbf X\mathbf v=\lVert\mathbf X\mathbf v\rVert^2$, so $\mathbf X\mathbf v=\mathbf 0$. Full column rank then forces $\mathbf v=\mathbf 0$. A square matrix with trivial null space is invertible (Ch. 8).

*Part 1 (minimiser).* Define $\hat{\boldsymbol\beta}$ by the formula. Then $\mathbf X^\mathsf T\mathbf X\hat{\boldsymbol\beta}=\mathbf X^\mathsf T\mathbf y$ (the **normal equations**), i.e. $\mathbf X^\mathsf T\mathbf e=\mathbf 0$. For any $\boldsymbol\beta$, write $\mathbf y-\mathbf X\boldsymbol\beta=\mathbf e+\mathbf X(\hat{\boldsymbol\beta}-\boldsymbol\beta)$ and expand:

$$
S(\boldsymbol\beta)=\lVert\mathbf e\rVert^2+2(\hat{\boldsymbol\beta}-\boldsymbol\beta)^\mathsf T\mathbf X^\mathsf T\mathbf e+\lVert\mathbf X(\hat{\boldsymbol\beta}-\boldsymbol\beta)\rVert^2=S(\hat{\boldsymbol\beta})+\lVert\mathbf X(\hat{\boldsymbol\beta}-\boldsymbol\beta)\rVert^2 .
$$

The last term is $\ge0$, with equality iff $\mathbf X(\hat{\boldsymbol\beta}-\boldsymbol\beta)=\mathbf 0$, iff $\boldsymbol\beta=\hat{\boldsymbol\beta}$ (full rank again). So $\hat{\boldsymbol\beta}$ is the unique minimiser. This argument needs no calculus. It is Pythagoras' theorem for the orthogonal projection $\hat{\mathbf y}=\mathbf X\hat{\boldsymbol\beta}$ of $\mathbf y$ onto the column space of $\mathbf X$ (Ch. 9).

*Part 1 (simple regression).* The first row of $\mathbf X^\mathsf T\mathbf e=\mathbf 0$ is $\sum e_i=0$, i.e. $\bar y=\hat\beta_0+\hat\beta_1\bar x$. The second row is $\sum x_ie_i=0$. Subtract $\bar x$ times the first row to get $\sum(x_i-\bar x)e_i=0$. Substituting $e_i=(y_i-\bar y)-\hat\beta_1(x_i-\bar x)$ gives $S_{xy}-\hat\beta_1S_{xx}=0$. The two columns $(1,\dots,1)^\mathsf T$ and $(x_1,\dots,x_n)^\mathsf T$ are linearly independent iff the $x_i$ are not all equal, iff $S_{xx}>0$.

*Part 2.* Let $\mathbf A=(\mathbf X^\mathsf T\mathbf X)^{-1}\mathbf X^\mathsf T$, so that $\mathbf A\mathbf X=\mathbf I_p$ and $\hat{\boldsymbol\beta}=\mathbf A\mathbf y=\boldsymbol\beta+\mathbf A\boldsymbol\varepsilon$. Linearity gives $E[\hat{\boldsymbol\beta}]=\boldsymbol\beta$. For any deterministic matrix $\mathbf M$, $\operatorname{Cov}(\mathbf M\boldsymbol\varepsilon)=\mathbf M\operatorname{Cov}(\boldsymbol\varepsilon)\mathbf M^\mathsf T$. This is E.2, part 3, applied to each linear combination. So

$$
\operatorname{Cov}(\hat{\boldsymbol\beta})=\sigma^2\mathbf A\mathbf A^\mathsf T=\sigma^2(\mathbf X^\mathsf T\mathbf X)^{-1}\mathbf X^\mathsf T\mathbf X(\mathbf X^\mathsf T\mathbf X)^{-1}=\sigma^2(\mathbf X^\mathsf T\mathbf X)^{-1}.
$$

For the residuals, let $\mathbf H=\mathbf X\mathbf A$ (the hat matrix). It is symmetric and idempotent ($\mathbf H^2=\mathbf H$), and $(\mathbf I-\mathbf H)\mathbf X=\mathbf 0$. Hence $\mathbf e=(\mathbf I-\mathbf H)\mathbf y=(\mathbf I-\mathbf H)\boldsymbol\varepsilon$ and $\mathrm{SSE}=\boldsymbol\varepsilon^\mathsf T(\mathbf I-\mathbf H)\boldsymbol\varepsilon$. For any square $\mathbf M$, $E[\boldsymbol\varepsilon^\mathsf T\mathbf M\boldsymbol\varepsilon]=\sum_{i,k}M_{ik}E[\varepsilon_i\varepsilon_k]=\sigma^2\operatorname{tr}\mathbf M$. Using $\operatorname{tr}(\mathbf{PQ})=\operatorname{tr}(\mathbf{QP})$,

$$
E[\mathrm{SSE}]=\sigma^2\operatorname{tr}(\mathbf I_n-\mathbf H)=\sigma^2\big(n-\operatorname{tr}\big((\mathbf X^\mathsf T\mathbf X)^{-1}\mathbf X^\mathsf T\mathbf X\big)\big)=\sigma^2(n-\operatorname{tr}\mathbf I_p)=(n-p)\sigma^2 .
$$

For simple regression, $\mathbf X^\mathsf T\mathbf X=\begin{pmatrix}n&\sum x_i\\\sum x_i&\sum x_i^2\end{pmatrix}$ has determinant $n\sum x_i^2-(\sum x_i)^2=nS_{xx}$. Its inverse is $\frac1{nS_{xx}}\begin{pmatrix}\sum x_i^2&-\sum x_i\\-\sum x_i&n\end{pmatrix}$, whose diagonal gives the two variances.

*Part 3 (Gauss–Markov).* Unbiasedness for every $\boldsymbol\beta$ means $E[\mathbf C\mathbf y]=\mathbf C\mathbf X\boldsymbol\beta=\boldsymbol\beta$ for all $\boldsymbol\beta$, i.e. $\mathbf C\mathbf X=\mathbf I_p$. Put $\mathbf D=\mathbf C-\mathbf A$. Then $\mathbf D\mathbf X=\mathbf I_p-\mathbf I_p=\mathbf 0$, and therefore

$$
\mathbf A\mathbf D^\mathsf T=(\mathbf X^\mathsf T\mathbf X)^{-1}(\mathbf D\mathbf X)^\mathsf T=\mathbf 0,\qquad\mathbf D\mathbf A^\mathsf T=(\mathbf A\mathbf D^\mathsf T)^\mathsf T=\mathbf 0 .
$$

Since $\tilde{\boldsymbol\beta}=\mathbf C\mathbf X\boldsymbol\beta+\mathbf C\boldsymbol\varepsilon$,

$$
\operatorname{Cov}(\tilde{\boldsymbol\beta})=\sigma^2\mathbf C\mathbf C^\mathsf T=\sigma^2(\mathbf A+\mathbf D)(\mathbf A+\mathbf D)^\mathsf T=\sigma^2\big(\mathbf A\mathbf A^\mathsf T+\mathbf D\mathbf D^\mathsf T\big)=\operatorname{Cov}(\hat{\boldsymbol\beta})+\sigma^2\mathbf D\mathbf D^\mathsf T .
$$

Finally, $\mathbf a^\mathsf T\mathbf D\mathbf D^\mathsf T\mathbf a=\lVert\mathbf D^\mathsf T\mathbf a\rVert^2\ge0$, so the difference is positive semidefinite, and $\operatorname{Var}(\mathbf a^\mathsf T\tilde{\boldsymbol\beta})-\operatorname{Var}(\mathbf a^\mathsf T\hat{\boldsymbol\beta})=\sigma^2\lVert\mathbf D^\mathsf T\mathbf a\rVert^2\ge0$. Equality for every $\mathbf a$ holds iff $\mathbf D=\mathbf 0$, so least squares is the *unique* BLUE. ∎

**Remarks**

- *No Gaussian assumption.* Gauss–Markov uses only the first two moments of $\boldsymbol\varepsilon$. With Gaussian errors, $\hat{\boldsymbol\beta}$ is also the maximum-likelihood estimator and attains the Cramér–Rao bound, so it is best among *all* unbiased estimators, linear or not ([Section 20.6](../part-6-probability-statistics/20-statistical-analysis.md)). Without Gaussianity, nonlinear estimators can win. For heavy-tailed (e.g. Laplacian) errors, least-absolute-deviation fitting has smaller variance than least squares.
- *Each hypothesis matters.*
  - **Equal variance.** If $\operatorname{Cov}(\boldsymbol\varepsilon)=\sigma^2\mathbf V$ with known $\mathbf V\ne\mathbf I$, then ordinary least squares is still unbiased but no longer best. Applying the theorem to the whitened model $\mathbf V^{-1/2}\mathbf y=\mathbf V^{-1/2}\mathbf X\boldsymbol\beta+\mathbf V^{-1/2}\boldsymbol\varepsilon$ shows that the BLUE is generalised least squares, $(\mathbf X^\mathsf T\mathbf V^{-1}\mathbf X)^{-1}\mathbf X^\mathsf T\mathbf V^{-1}\mathbf y$. With diagonal $\mathbf V$ this is the weighted least squares of [Section 20.9](../part-6-probability-statistics/20-statistical-analysis.md).
  - **Unbiasedness.** Biased linear estimators can have smaller mean-square error. Ridge regression $(\mathbf X^\mathsf T\mathbf X+\lambda\mathbf I)^{-1}\mathbf X^\mathsf T\mathbf y$ does, for suitable $\lambda>0$, when $\mathbf X^\mathsf T\mathbf X$ is ill-conditioned.
  - **Full rank.** With collinear regressors, $\boldsymbol\beta$ is not identifiable. Many $\boldsymbol\beta$ give the same $\mathbf X\boldsymbol\beta$ ([Section 25.8](../part-7-numerical-computational/25-mathematical-modeling.md)).
- *Divisor $n-p$.* Part 2 shows that $s^2=\mathrm{SSE}/(n-p)$ is unbiased for $\sigma^2$. E.6 is the special case $p=1$, $\mathbf X=(1,\dots,1)^\mathsf T$, $\hat\beta=\bar y$.
- *Numerics.* The formula $(\mathbf X^\mathsf T\mathbf X)^{-1}\mathbf X^\mathsf T\mathbf y$ is for analysis. Compute $\hat{\boldsymbol\beta}$ by QR factorisation ([Section 24.3](../part-7-numerical-computational/24-computational-mathematics.md)), because $\kappa(\mathbf X^\mathsf T\mathbf X)=\kappa(\mathbf X)^2$.
- *Where used.* [Section 20.9](../part-6-probability-statistics/20-statistical-analysis.md), discrete least squares in [Section 22.6](../part-7-numerical-computational/22-numerical-methods.md), and parameter identification in [Section 25.8](../part-7-numerical-computational/25-mathematical-modeling.md).

## E.8 The law of propagation of uncertainty (GUM)

> [!definition] Measurement model
> A measurand is modelled as $Y=f(\mathbf X)$, where $\mathbf X=(X_1,\dots,X_N)^\mathsf T$ is a random vector of input quantities. The input **estimates** are $\mathbf x=E[\mathbf X]$, the **standard uncertainties** are $u(x_i)=\sqrt{\operatorname{Var}X_i}$, and the **covariances** are $u(x_i,x_k)=\operatorname{Cov}(X_i,X_k)$. Together they form the input covariance matrix $\boldsymbol\Sigma_x$. The **sensitivity coefficients** are $c_i=\partial f/\partial x_i$ evaluated at $\mathbf x$, and $\mathbf c=\nabla f(\mathbf x)$.

**Statement**

> [!theorem] First-order propagation of uncertainty and its accuracy
> Let $f:\mathbb R^N\to\mathbb R$ be twice continuously differentiable, with Hessian bounded in the spectral norm, $\lVert\nabla^2f(\boldsymbol\xi)\rVert_2\le M$ for all $\boldsymbol\xi$. Let $\mathbf X$ have mean $\mathbf x$, covariance $\boldsymbol\Sigma_x$, and finite fourth moments, $m_4=E\lVert\mathbf X-\mathbf x\rVert^4<\infty$. Define the linearised model $Y_{\rm lin}=f(\mathbf x)+\mathbf c^\mathsf T(\mathbf X-\mathbf x)$ and the **combined standard uncertainty**
> $$
> u_c^2(y)=\mathbf c^\mathsf T\boldsymbol\Sigma_x\mathbf c=\sum_{i=1}^Nc_i^2u^2(x_i)+2\sum_{i=1}^{N-1}\sum_{k=i+1}^Nc_ic_k\,u(x_i,x_k).
> $$
> Then:
>
> 1. $E[Y_{\rm lin}]=f(\mathbf x)$ and $\operatorname{Var}(Y_{\rm lin})=u_c^2(y)$ **exactly**;
> 2. $\lvert Y-Y_{\rm lin}\rvert\le\frac M2\lVert\mathbf X-\mathbf x\rVert^2$ at every outcome, hence $\lvert E[Y]-f(\mathbf x)\rvert\le\frac M2\operatorname{tr}\boldsymbol\Sigma_x$;
> 3. the true standard deviation $\sigma_Y$ of $Y$ satisfies
> $$
> \big\lvert\sigma_Y-u_c(y)\big\rvert\le\frac M2\sqrt{m_4}.
> $$
>
> If all input deviations scale with a small parameter $u$ ($\mathbf X-\mathbf x=u\mathbf W$ with $\mathbf W$ fixed), then $u_c=O(u)$ while the errors in parts 2 and 3 are $O(u^2)$. The relative error of the GUM formula is therefore $O(u)$, provided $\mathbf c\ne\mathbf 0$.

**Proof**

*Part 1.* $Y_{\rm lin}$ is the constant $f(\mathbf x)-\mathbf c^\mathsf T\mathbf x$ plus the linear form $\mathbf c^\mathsf T\mathbf X$. By linearity of expectation, $E[Y_{\rm lin}]=f(\mathbf x)+\mathbf c^\mathsf T(E\mathbf X-\mathbf x)=f(\mathbf x)$. Adding a constant does not change a variance, so E.2, part 3, gives $\operatorname{Var}(Y_{\rm lin})=\mathbf c^\mathsf T\boldsymbol\Sigma_x\mathbf c$. Written out, the double sum $\sum_{i,k}c_ic_k(\boldsymbol\Sigma_x)_{ik}$ has diagonal terms $c_i^2u^2(x_i)$ and each off-diagonal pair $(i,k)$, $(k,i)$ appearing twice. That is the displayed formula.

*Part 2 (Taylor's theorem).* Fix an outcome and put $\mathbf d=\mathbf X-\mathbf x$ and $g(t)=f(\mathbf x+t\mathbf d)$ for $t\in[0,1]$. Then $g\in C^2$, $g'(t)=\nabla f(\mathbf x+t\mathbf d)^\mathsf T\mathbf d$ and $g''(t)=\mathbf d^\mathsf T\nabla^2f(\mathbf x+t\mathbf d)\,\mathbf d$ (chain rule, Ch. 6). Taylor's theorem with Lagrange remainder (Ch. 4) gives $g(1)=g(0)+g'(0)+\frac12g''(\theta)$ for some $\theta\in(0,1)$, i.e.

$$
R:=Y-Y_{\rm lin}=\tfrac12\,\mathbf d^\mathsf T\nabla^2f(\mathbf x+\theta\mathbf d)\,\mathbf d,\qquad\lvert R\rvert\le\tfrac12\lVert\nabla^2f\rVert_2\lVert\mathbf d\rVert^2\le\tfrac M2\lVert\mathbf X-\mathbf x\rVert^2 .
$$

The quadratic-form bound $\lvert\mathbf d^\mathsf T\mathbf B\mathbf d\rvert\le\lVert\mathbf B\rVert_2\lVert\mathbf d\rVert^2$ is the Cauchy–Schwarz inequality (Ch. 9). Taking expectations, $\lvert E[Y]-f(\mathbf x)\rvert=\lvert E[R]\rvert\le E\lvert R\rvert\le\frac M2E\lVert\mathbf X-\mathbf x\rVert^2=\frac M2\operatorname{tr}\boldsymbol\Sigma_x$.

*Part 3.* First, $Y$ has a finite second moment, because $\lvert Y\rvert\le\lvert f(\mathbf x)\rvert+\lVert\mathbf c\rVert\lVert\mathbf d\rVert+\frac M2\lVert\mathbf d\rVert^2$ and $m_4<\infty$. Write $\lVert Z\rVert_{2}=\sqrt{E[Z^2]}$, which is a norm on random variables with finite second moment (Minkowski's inequality). Then $\sigma_Y=\lVert Y-EY\rVert_2$ and $u_c=\lVert Y_{\rm lin}-f(\mathbf x)\rVert_2$, and

$$
Y-E[Y]=\big(Y_{\rm lin}-f(\mathbf x)\big)+\big(R-E[R]\big).
$$

By the reverse triangle inequality for the norm,

$$
\lvert\sigma_Y-u_c\rvert\le\lVert R-E[R]\rVert_2=\sqrt{\operatorname{Var}R}\le\sqrt{E[R^2]}\le\frac M2\sqrt{E\lVert\mathbf X-\mathbf x\rVert^4}.
$$

*Scaling.* With $\mathbf d=u\mathbf W$, $u_c=u\sqrt{\mathbf c^\mathsf T\operatorname{Cov}(\mathbf W)\mathbf c}$ is proportional to $u$, while $\operatorname{tr}\boldsymbol\Sigma_x$ and $\sqrt{m_4}$ are proportional to $u^2$. ∎

**Remarks**

- *Products and quotients.* For $Y=A\,X_1^{p_1}\cdots X_N^{p_N}$ with $x_i\ne0$, $c_i=p_iy/x_i$. For uncorrelated inputs, part 1 becomes the relative-uncertainty rule $(u_c/y)^2=\sum(p_iu(x_i)/x_i)^2$ of [Section 20.10](../part-6-probability-statistics/20-statistical-analysis.md). Such $f$ do not have a globally bounded Hessian. The theorem then applies with $M$ the bound over a region that contains $\mathbf X$ with overwhelming probability, and the part of the distribution outside that region must be handled separately. This is the usual situation in practice.
- *Counterexample: vanishing sensitivity.* Let $Y=X^2$ with $X\sim\mathcal N(x,u^2)$. Then exactly $\operatorname{Var}Y=4x^2u^2+2u^4$, while $u_c=2\lvert x\rvert u$. At $x=1$, $u=0.1$: $\sigma_Y=0.20050$ against $u_c=0.2$, well inside the bound $\frac M2\sqrt{m_4}=\sqrt3\,u^2=0.0173$ ($M=2$, $m_4=3u^4$). At $x=0$, however, $c=0$ and $u_c=0$, yet $\sigma_Y=\sqrt2\,u^2$. The first-order term has vanished, and the "$O(u^2)$ error" is the whole answer. The GUM then requires second-order terms or Monte Carlo propagation ([Section 24.8](../part-7-numerical-computational/24-computational-mathematics.md)).
- *Second-order mean correction.* The proof also shows that $E[Y]=f(\mathbf x)+\frac12\operatorname{tr}\big(\nabla^2f(\mathbf x)\boldsymbol\Sigma_x\big)+o(u^2)$ when $\nabla^2f$ is continuous. A curved model *biases* the output. For example, $E[X^2]=x^2+u^2$ exactly.
- *Correlation.* Part 1 is an identity for the linearised model, so the covariance terms are not optional. With $r=\pm1$ they turn root-sum-of-squares addition into linear addition or cancellation, as the chapter discusses.
- *Where used.* [Section 20.10](../part-6-probability-statistics/20-statistical-analysis.md), propagation of errors in [Section 22.2](../part-7-numerical-computational/22-numerical-methods.md), and linearised uncertainty propagation in [Section 24.8](../part-7-numerical-computational/24-computational-mathematics.md).

## E.9 The Wiener–Khinchin theorem and filtered noise

> [!definition] Power spectral density
> Let $X(t)$ be a real wide-sense stationary (WSS) process with autocorrelation $R_{XX}(\tau)=E[X(t+\tau)X(t)]$. Assume $X$ is **mean-square continuous** ($R_{XX}$ is continuous at $0$, hence everywhere, Section 21.3) and has jointly measurable sample paths. For $T>0$, let $X_T(\omega)=\int_{-T}^{T}X(t)e^{-j\omega t}\,dt$. The **power spectral density** is
> $$
> S_{XX}(\omega)=\lim_{T\to\infty}\frac{E\big[\lvert X_T(\omega)\rvert^2\big]}{2T},
> $$
> whenever the limit exists.

**Statement**

> [!theorem] Wiener–Khinchin
> Under the assumptions of the definition, suppose $R_{XX}\in L^1(\mathbb R)$, i.e. $\int_{-\infty}^\infty\lvert R_{XX}(\tau)\rvert\,d\tau<\infty$. Then for every $\omega$ the limit defining $S_{XX}(\omega)$ exists and
> $$
> S_{XX}(\omega)=\int_{-\infty}^{\infty}R_{XX}(\tau)\,e^{-j\omega\tau}\,d\tau .
> $$
> $S_{XX}$ is real, even, non-negative, continuous and bounded by $\lVert R_{XX}\rVert_1$. If also $S_{XX}\in L^1$, then for every $\tau$
> $$
> R_{XX}(\tau)=\frac1{2\pi}\int_{-\infty}^\infty S_{XX}(\omega)e^{j\omega\tau}\,d\omega,\qquad E[X^2(t)]=\frac1{2\pi}\int_{-\infty}^\infty S_{XX}(\omega)\,d\omega .
> $$

> [!theorem] Lemma E.9a (diagonal integration)
> If $g$ is integrable on $[-2T,2T]$, then
> $$
> \int_{-T}^{T}\!\!\int_{-T}^{T}g(t-s)\,dt\,ds=\int_{-2T}^{2T}\big(2T-\lvert\tau\rvert\big)g(\tau)\,d\tau .
> $$

*Proof.* In the inner integral substitute $\tau=t-s$ for fixed $s$. This gives $\int_{-T}^{T}\int_{-T-s}^{T-s}g(\tau)\,d\tau\,ds$, the integral of $g(\tau)$ over the region $\{(\tau,s):\lvert s\rvert\le T,\ \lvert s+\tau\rvert\le T\}$. By Fubini's theorem, which applies because $\iint\lvert g(t-s)\rvert\,dt\,ds\le2T\int_{-2T}^{2T}\lvert g\rvert<\infty$, we may integrate over $s$ first. For fixed $\tau$ with $\lvert\tau\rvert\le2T$, the admissible $s$ form the interval $[\max(-T,-T-\tau),\ \min(T,T-\tau)]$. Its length is $2T-\lvert\tau\rvert$. ∎

**Proof of the theorem**

*Step 1 ($X_T$ exists and $E\lvert X_T\rvert^2$ is a double integral).* By Cauchy–Schwarz, $E\lvert X(t)X(s)\rvert\le\sqrt{E[X^2(t)]E[X^2(s)]}=R_{XX}(0)$. Therefore

$$
E\int_{-T}^{T}\!\!\int_{-T}^{T}\lvert X(t)X(s)\rvert\,dt\,ds\le4T^2R_{XX}(0)<\infty .
$$

So, with probability one, $X$ is integrable on $[-T,T]$, $X_T(\omega)$ exists, and Fubini's theorem allows expectation and integration to be interchanged:

$$
E\lvert X_T(\omega)\rvert^2=E\Big[\int_{-T}^{T}\!\!\int_{-T}^{T}X(t)X(s)e^{-j\omega(t-s)}\,dt\,ds\Big]=\int_{-T}^{T}\!\!\int_{-T}^{T}R_{XX}(t-s)\,e^{-j\omega(t-s)}\,dt\,ds .
$$

Here $\lvert X_T\rvert^2=X_T\overline{X_T}$ was used, and $X$ is real.

*Step 2 (reduce to one integral).* Lemma E.9a with $g(\tau)=R_{XX}(\tau)e^{-j\omega\tau}$ gives

$$
\frac{E\lvert X_T(\omega)\rvert^2}{2T}=\int_{-\infty}^{\infty}w_T(\tau)\,R_{XX}(\tau)e^{-j\omega\tau}\,d\tau,\qquad w_T(\tau)=\Big(1-\frac{\lvert\tau\rvert}{2T}\Big)\mathbf 1_{\{\lvert\tau\rvert\le2T\}} .
$$

*Step 3 (limit).* For every $\tau$, $w_T(\tau)\to1$ as $T\to\infty$, and $\lvert w_TR_{XX}e^{-j\omega\tau}\rvert\le\lvert R_{XX}(\tau)\rvert$, which is integrable. The dominated convergence theorem gives the stated limit.

*Step 4 (properties).* $R_{XX}$ is even, since $R_{XX}(-\tau)=E[X(t-\tau)X(t)]=R_{XX}(\tau)$ by stationarity. So $S_{XX}(\omega)=\int R_{XX}(\tau)\cos\omega\tau\,d\tau$ is real and even. It is non-negative as a limit of the non-negative numbers $E\lvert X_T\rvert^2/(2T)$. Also $\lvert S_{XX}\rvert\le\lVert R_{XX}\rVert_1$, and continuity in $\omega$ follows by dominated convergence again.

*Step 5 (inversion).* If $S_{XX}\in L^1$, the Fourier inversion theorem for $L^1$ functions with $L^1$ transforms (Ch. 16) gives $R_{XX}(\tau)=\frac1{2\pi}\int S_{XX}e^{j\omega\tau}d\omega$ for almost every $\tau$. Both sides are continuous in $\tau$ (the right by dominated convergence), so equality holds for every $\tau$. At $\tau=0$ we get $E[X^2]=R_{XX}(0)$. ∎

```desmos-graph
left=-12; right=12; bottom=-0.05; top=1.05
---
y=e^{-\operatorname{abs}\left(x\right)}|black|dashed
y=\left(1-\frac{\operatorname{abs}\left(x\right)}{2}\right)e^{-\operatorname{abs}\left(x\right)}|red|-2<x<2
y=\left(1-\frac{\operatorname{abs}\left(x\right)}{10}\right)e^{-\operatorname{abs}\left(x\right)}|blue|-10<x<10
y=1-\frac{\operatorname{abs}\left(x\right)}{2}|red|dotted|-2<x<2
y=1-\frac{\operatorname{abs}\left(x\right)}{10}|blue|dotted|-10<x<10
```

*Step 2 in pictures for $R_{XX}(\tau)=e^{-\lvert\tau\rvert}$ (dashed): a record of length $2T$ sees the autocorrelation through the triangular window $w_T$ (dotted), giving the products shown for $T=1$ (red) and $T=5$ (blue). As $T\to\infty$ the window flattens to $1$. At $\omega=0.5$, $E\lvert X_T\rvert^2/2T=1.082,\ 1.504,\ 1.590$ for $T=1,5,50$, against $S_{XX}(0.5)=2/(1+0.25)=1.6$.*

> [!theorem] Random signals through LTI systems
> Let $X$ satisfy the hypotheses of the Wiener–Khinchin theorem. Let $h$ be real with $h\in L^1(\mathbb R)$, i.e. the system is BIBO stable (Ch. 14), and let $H(\omega)=\int h(t)e^{-j\omega t}\,dt$. Define $Y(t)=\int_{-\infty}^{\infty}h(\alpha)X(t-\alpha)\,d\alpha$. Then:
>
> 1. $Y$ is a well-defined WSS process with $E[Y]=\mu_XH(0)$ and
> $$
> R_{YY}(\tau)=\int_{-\infty}^{\infty}\!\!\int_{-\infty}^{\infty}h(\alpha)h(\beta)\,R_{XX}(\tau-\alpha+\beta)\,d\alpha\,d\beta ,
> $$
> a continuous function in $L^1$ with $\lVert R_{YY}\rVert_1\le\lVert h\rVert_1^2\lVert R_{XX}\rVert_1$;
> 2. $S_{YY}(\omega)=\lvert H(\omega)\rvert^2S_{XX}(\omega)$ for every $\omega$;
> 3. if $S_{XX}\in L^1$, then $E[Y^2]=\dfrac1{2\pi}\displaystyle\int_{-\infty}^{\infty}\lvert H(\omega)\rvert^2S_{XX}(\omega)\,d\omega$.

**Proof**

*Part 1.* By Cauchy–Schwarz, $E\lvert X(t)\rvert\le\sqrt{R_{XX}(0)}$, so $E\int\lvert h(\alpha)\rvert\lvert X(t-\alpha)\rvert\,d\alpha\le\lVert h\rVert_1\sqrt{R_{XX}(0)}<\infty$. Hence the integral defining $Y(t)$ converges absolutely with probability one. By Fubini, $E[Y(t)]=\int h(\alpha)\mu_X\,d\alpha=\mu_XH(0)$. For the correlation, the integrand of

$$
E[Y(t+\tau)Y(t)]=E\Big[\iint h(\alpha)h(\beta)X(t+\tau-\alpha)X(t-\beta)\,d\alpha\,d\beta\Big]
$$

is dominated in expectation by $\lvert h(\alpha)h(\beta)\rvert R_{XX}(0)$, whose double integral is $\lVert h\rVert_1^2R_{XX}(0)<\infty$. So Fubini again moves $E$ inside, and $E[X(t+\tau-\alpha)X(t-\beta)]=R_{XX}(\tau-\alpha+\beta)$. The result depends on $\tau$ only, and $E[Y^2]=R_{YY}(0)\le\lVert h\rVert_1^2R_{XX}(0)<\infty$, so $Y$ is WSS. Integrating $\lvert R_{YY}\rvert$ over $\tau$ and using Tonelli's theorem with the substitution $u=\tau-\alpha+\beta$ gives $\lVert R_{YY}\rVert_1\le\lVert h\rVert_1^2\lVert R_{XX}\rVert_1$. Continuity of $R_{YY}$ follows from continuity of $R_{XX}$, the bound $\lvert R_{XX}\rvert\le R_{XX}(0)$, and dominated convergence.

*Part 2.* By part 1, $Y$ satisfies the hypotheses of the Wiener–Khinchin theorem, so $S_{YY}(\omega)=\int R_{YY}(\tau)e^{-j\omega\tau}\,d\tau$. Insert the double integral. The triple integral converges absolutely (the bound in part 1), so Fubini permits any order of integration. Integrate over $\tau$ first with $u=\tau-\alpha+\beta$, noting $e^{-j\omega\tau}=e^{-j\omega u}e^{-j\omega\alpha}e^{j\omega\beta}$:

$$
S_{YY}(\omega)=\int h(\alpha)e^{-j\omega\alpha}d\alpha\cdot\int h(\beta)e^{j\omega\beta}d\beta\cdot\int R_{XX}(u)e^{-j\omega u}du=H(\omega)\,\overline{H(\omega)}\,S_{XX}(\omega).
$$

The middle factor is $\overline{H(\omega)}$ because $h$ is real. Hence $S_{YY}=\lvert H\rvert^2S_{XX}$.

*Part 3.* $\lvert H(\omega)\rvert\le\lVert h\rVert_1$, so $S_{YY}\le\lVert h\rVert_1^2S_{XX}\in L^1$. The inversion part of the Wiener–Khinchin theorem at $\tau=0$ applies to $Y$. ∎

**Remarks**

- *Why $R_{XX}\in L^1$.* The hypothesis excludes components with a discrete spectrum. For the random-phase sinusoid $X(t)=A\cos(\omega_0t+\Theta)$, $\Theta\sim\mathcal U(0,2\pi)$, we have $R_{XX}=\frac{A^2}2\cos\omega_0\tau\notin L^1$. Then $E\lvert X_T(\omega_0)\rvert^2/(2T)$ grows like $T$, and the PSD exists only as the generalised function $\frac{\pi A^2}{2}[\delta(\omega-\omega_0)+\delta(\omega+\omega_0)]$. A non-zero mean behaves the same way, contributing $2\pi\mu_X^2\delta(\omega)$. The chapter's table of transform pairs ([Section 21.4](../part-6-probability-statistics/21-random-processes.md)) uses this distributional extension. The chapter's sufficient condition $\int\lvert\tau R_{XX}\rvert<\infty$ (together with $R_{XX}\in L^1$) gives the stronger conclusion that the convergence is $O(1/T)$, because $\big\lvert E\lvert X_T\rvert^2/(2T)-S_{XX}\big\rvert\le\frac1{2T}\int\lvert\tau R_{XX}\rvert\,d\tau+\int_{\lvert\tau\rvert>2T}\lvert R_{XX}\rvert\,d\tau$.
- *Non-negativity for free.* Step 4 shows that the Fourier transform of a valid autocorrelation is non-negative. This is the easy half of Bochner's theorem. For example, a rectangular "autocorrelation" $\mathbf 1_{\{\lvert\tau\rvert\le1\}}$ has transform $2\sin\omega/\omega$, which takes negative values, so no WSS process has it.
- *White noise.* $R=\frac{N_0}2\delta$ is not a function, so white noise lies outside the theorem. Its role is as a limit: band-limited white noise satisfies the hypotheses, and letting the bandwidth grow in part 2 gives $S_{YY}=\frac{N_0}{2}\lvert H\rvert^2$ for every stable $H$ with $\int\lvert H\rvert^2<\infty$. For the RC low-pass $H=1/(1+j\omega RC)$, part 3 gives $E[Y^2]=\frac1{2\pi}\int\frac{N_0/2}{1+\omega^2R^2C^2}\,d\omega=\frac{N_0}{4RC}$, confirmed by numerical quadrature. This matches the chapter's worked example.
- *Where used.* [Section 21.4](../part-6-probability-statistics/21-random-processes.md) (PSD, LTI systems with random inputs), noise-equivalent bandwidth and coloured noise in [Section 21.6](../part-6-probability-statistics/21-random-processes.md), and the matched filter in [Section 21.8](../part-6-probability-statistics/21-random-processes.md).

## E.10 Quadratic convergence of Newton–Raphson

Newton's iteration for $f(x)=0$ is $x_{k+1}=x_k-f(x_k)/f'(x_k)$ ([Section 22.4](../part-7-numerical-computational/22-numerical-methods.md)). Write $e_k=x_k-x^\ast$ for the error.

**Statement**

> [!theorem] Local quadratic convergence of Newton's method
> Let $f(x^\ast)=0$, and let $f$ be twice continuously differentiable on $I=[x^\ast-r,\,x^\ast+r]$ with $f'\ne0$ on $I$ (a **simple root**; if only $f'(x^\ast)\ne0$ is known, shrink $r$ using continuity of $f'$). Put
> $$
> m=\min_I\lvert f'\rvert>0,\qquad M=\max_I\lvert f''\rvert,\qquad K=\frac{M}{2m},\qquad\delta=\min\Big(r,\ \frac1{2K}\Big)\quad(\delta=r\text{ if }M=0).
> $$
> If $\lvert x_0-x^\ast\rvert\le\delta$, then all iterates are defined and lie in $I$, and
> 1. $\lvert e_{k+1}\rvert\le K\,e_k^2$ and $\lvert e_{k+1}\rvert\le\frac12\lvert e_k\rvert$ for every $k$, so $x_k\to x^\ast$;
> 2. $K\lvert e_k\rvert\le\big(K\lvert e_0\rvert\big)^{2^k}\le2^{-2^k}$: the number of correct digits roughly doubles at every step;
> 3. if $e_k\ne0$ for all $k$, then $\displaystyle\lim_{k\to\infty}\frac{e_{k+1}}{e_k^2}=\frac{f''(x^\ast)}{2f'(x^\ast)}$.

**Proof**

*Step 1 (one-step error identity).* Suppose $x_k\in I$. Then $f'(x_k)\ne0$, so $x_{k+1}$ is defined. Taylor's theorem with Lagrange remainder (Ch. 4), expanding about $x_k$ and evaluating at $x^\ast$, gives

$$
0=f(x^\ast)=f(x_k)+f'(x_k)(x^\ast-x_k)+\tfrac12f''(\xi_k)(x^\ast-x_k)^2
$$

for some $\xi_k$ between $x_k$ and $x^\ast$, hence $\xi_k\in I$. Divide by $f'(x_k)$ and use $f(x_k)/f'(x_k)=x_k-x_{k+1}$:

$$
0=(x_k-x_{k+1})+(x^\ast-x_k)+\frac{f''(\xi_k)}{2f'(x_k)}e_k^2\quad\Longrightarrow\quad e_{k+1}=\frac{f''(\xi_k)}{2f'(x_k)}\,e_k^2 .
$$

So $\lvert e_{k+1}\rvert\le K e_k^2$.

*Step 2 (induction: the iterates stay in $I$).* Suppose $\lvert e_k\rvert\le\delta$, which holds for $k=0$. Then $x_k\in I$, and by Step 1, $\lvert e_{k+1}\rvert\le(K\lvert e_k\rvert)\lvert e_k\rvert\le(K\delta)\lvert e_k\rvert\le\frac12\lvert e_k\rvert\le\delta$. By induction every iterate satisfies $\lvert e_k\rvert\le\delta$, and $\lvert e_k\rvert\le2^{-k}\lvert e_0\rvert\to0$. This proves part 1. (If some $e_k=0$, all later errors are $0$, since $f(x^\ast)=0$.)

*Step 3 (doubling of digits).* Let $q_k=K\lvert e_k\rvert$. Step 1 gives $q_{k+1}\le q_k^2$, so by induction $q_k\le q_0^{2^k}$, and $q_0\le K\delta\le\frac12$.

*Step 4 (the limit).* As $k\to\infty$, $x_k\to x^\ast$ and $\xi_k\to x^\ast$ (it lies between them). By continuity of $f''$ and $f'$, the factor $f''(\xi_k)/(2f'(x_k))$ tends to $f''(x^\ast)/(2f'(x^\ast))$. ∎

**Remarks**

- *Numerical check.* For $f(x)=x^3-x-2$ (Worked Example 22.14), the root is $x^\ast=1.5213797$. The predicted constant is $f''(x^\ast)/(2f'(x^\ast))=6x^\ast/(2(3x^{\ast2}-1))=0.768$, matching the observed ratios $0.79,\ 0.77,\ 0.77$ in the chapter's table.
- *Counterexample 1: multiple roots.* For $f(x)=x^2$ (double root, $f'(0)=0$) Newton gives $x_{k+1}=x_k/2$. The convergence is only linear with rate $\frac12$. In general, a root of multiplicity $p$ gives the linear rate $1-1/p$ ([Section 22.4](../part-7-numerical-computational/22-numerical-methods.md), "Multiple roots"). The hypothesis $f'(x^\ast)\ne0$ is what makes $m>0$.
- *Counterexample 2: the result is local.* For $f(x)=\arctan x$, the iteration is $x_{k+1}=x_k-(1+x_k^2)\arctan x_k$. It converges for $\lvert x_0\rvert<x_c$ and diverges with growing oscillation for $\lvert x_0\rvert>x_c$. Here $x_c=1.39175$ solves $2x=(1+x^2)\arctan x$. From $x_0=1.5$ the iterates are $-1.694,\ 2.321,\ -5.114,\ 32.3,\dots$. From $x_0=1.3$ they are $-1.162,\ 0.859,\ -0.374,\ 0.034,\ -2.6\times10^{-5},\ 1.2\times10^{-14}$. The theorem's $\delta$ is a safe radius, not the exact basin.
- *Smoothness.* Suppose that instead of $f\in C^2$ we only know that $f'$ is Hölder continuous with exponent $\alpha\in(0,1]$. Writing the Taylor remainder as $\int_{x_k}^{x^\ast}(f'(s)-f'(x_k))\,ds=O(\lvert e_k\rvert^{1+\alpha})$, Step 1 gives order $1+\alpha$ instead of $2$. At the other extreme, if $f''(x^\ast)=0$ and $f\in C^3$, the same expansion carried one order further gives cubic convergence.
- *Where used.* [Section 22.4](../part-7-numerical-computational/22-numerical-methods.md) (including diode circuit analysis), Gauss–Newton in [Section 24.7](../part-7-numerical-computational/24-computational-mathematics.md), and Newton steps inside implicit ODE solvers ([Section 22.9](../part-7-numerical-computational/22-numerical-methods.md)).

## E.11 The bisection error bound

Bisection starts from $[a_0,b_0]$ with $f(a_0)f(b_0)<0$. At step $k$ it sets $m_k=\frac12(a_k+b_k)$. It stops if $f(m_k)=0$. Otherwise it keeps the half on which $f$ changes sign ([Section 22.4](../part-7-numerical-computational/22-numerical-methods.md)).

**Statement**

> [!theorem] Convergence of bisection
> Let $f$ be continuous on $[a_0,b_0]$ with $f(a_0)f(b_0)<0$. If the algorithm never stops at an exact zero, then:
> 1. $a_0\le a_1\le\cdots\le a_k<b_k\le\cdots\le b_1\le b_0$, $f(a_k)f(b_k)<0$, and $b_k-a_k=(b_0-a_0)/2^k$;
> 2. $a_k$ and $b_k$ converge to a common limit $x^\ast$ with $f(x^\ast)=0$, and $a_k\le x^\ast\le b_k$ for every $k$;
> 3. $\displaystyle\lvert m_k-x^\ast\rvert\le\frac{b_0-a_0}{2^{k+1}}$ for every $k$. Hence $\lvert m_k-x^\ast\rvert\le\tau$ is guaranteed as soon as $k\ge\log_2\dfrac{b_0-a_0}{\tau}-1$.

**Proof**

*Part 1 (induction).* The case $k=0$ is the hypothesis. Suppose part 1 holds at step $k$. Since $f(m_k)\ne0$ and $f(a_k)f(b_k)<0$, exactly one of $f(a_k)f(m_k)<0$ and $f(m_k)f(b_k)<0$ holds: $f(m_k)$ has the sign of one endpoint value and the opposite sign of the other. The algorithm keeps that half. So $[a_{k+1},b_{k+1}]$ is either $[a_k,m_k]$ or $[m_k,b_k]$. In both cases the sign condition holds, the new interval is nested in the old one, and its width is $\frac12(b_k-a_k)=(b_0-a_0)/2^{k+1}$.

*Part 2.* $(a_k)$ is non-decreasing and bounded above by $b_0$. $(b_k)$ is non-increasing and bounded below by $a_0$. Monotone bounded sequences converge (completeness of $\mathbb R$, Ch. 4). Since $b_k-a_k\to0$, the limits coincide; call the common limit $x^\ast$. Monotonicity gives $a_k\le x^\ast\le b_k$. By continuity, $f(a_k)f(b_k)\to f(x^\ast)^2$. Each term is negative, so the limit satisfies $f(x^\ast)^2\le0$, i.e. $f(x^\ast)=0$.

*Part 3.* $x^\ast$ and $m_k$ both lie in $[a_k,b_k]$, and $m_k$ is its centre. So $\lvert m_k-x^\ast\rvert\le\frac12(b_k-a_k)=(b_0-a_0)/2^{k+1}$. The iteration count follows by solving $(b_0-a_0)/2^{k+1}\le\tau$ for $k$. ∎

**Remarks**

- *Stronger than "some root".* Part 2 identifies one root $x^\ast$ that is the limit and that satisfies the bound at *every* step. The intermediate value theorem by itself would only give a root in each bracket, possibly a different one each time.
- *Counterexample: continuity is essential.* $f(x)=1/x$ on $[-1,2]$ has $f(-1)f(2)<0$ but no root. Bisection keeps the sign change and converges to the pole $x=0$ (brackets $[-1,0.5]$, $[-0.25,0.5]$, $[-0.25,0.125]$, …). Part 2 fails exactly at the step "by continuity".
- *Sharpness.* The bound cannot be improved in general. If the root sits arbitrarily close to an endpoint of the current bracket, the midpoint is almost half a width away.
- *Floating point.* In finite precision the bracket stops shrinking at about one unit in the last place of $x^\ast$. Correct sign evaluation is all that is needed for the bound to hold until then ([Section 22.1](../part-7-numerical-computational/22-numerical-methods.md)).
- *Where used.* [Section 22.4](../part-7-numerical-computational/22-numerical-methods.md) (Example 22.13). It is also the safeguard inside hybrid root finders that fall back to bisection when Newton or secant steps leave the bracket.

## E.12 Error terms of the trapezoidal and Simpson rules

Throughout, the error of a rule $Q$ for $I=\int_a^bf\,dx$ is defined as $E=I-Q$.

**Statement**

> [!theorem] Trapezoidal and Simpson error terms
> 1. **Trapezoidal rule.** If $f\in C^2[a,b]$, then for some $\xi\in(a,b)$
> $$
> \int_a^bf\,dx-\frac{b-a}2\big(f(a)+f(b)\big)=-\frac{(b-a)^3}{12}f''(\xi).
> $$
> 2. **Simpson's rule.** If $f\in C^4[a,b]$, $h=\frac{b-a}2$ and $m=\frac{a+b}2$, then for some $\xi\in(a,b)$
> $$
> \int_a^bf\,dx-\frac h3\big(f(a)+4f(m)+f(b)\big)=-\frac{h^5}{90}f^{(4)}(\xi).
> $$
> 3. **Composite rules.** With $h=(b-a)/n$ and nodes $x_i=a+ih$, for some $\xi,\xi'\in(a,b)$
> $$
> I-T_n=-\frac{(b-a)h^2}{12}f''(\xi)\quad(f\in C^2),\qquad I-S_n=-\frac{(b-a)h^4}{180}f^{(4)}(\xi')\quad(f\in C^4,\ n\text{ even}).
> $$

Two lemmas carry the proof.

> [!theorem] Lemma E.12a (weighted mean value theorem)
> Let $\phi$ be continuous on $[a,b]$, and let $w\le0$ be integrable on $[a,b]$. Let $q$ be a function with $q\,w$ integrable and $\min\phi\le q(x)\le\max\phi$ wherever $w(x)\ne0$. Then $\int_a^bq\,w\,dx=\phi(\xi)\int_a^bw\,dx$ for some $\xi\in[a,b]$.

*Proof.* Let $\underline\phi=\min\phi$ and $\overline\phi=\max\phi$. Since $w\le0$, we have $\overline\phi\,w\le q\,w\le\underline\phi\,w$ pointwise. Integrating gives $\overline\phi W\le\int qw\le\underline\phi W$ with $W=\int w\le0$. If $W=0$, then $w=0$ almost everywhere, both sides vanish, and any $\xi$ works. Otherwise divide by $W<0$, which reverses the inequalities: $\underline\phi\le\int qw/W\le\overline\phi$. The intermediate value theorem then gives $\xi$ with $\phi(\xi)=\int qw/W$. ∎

> [!theorem] Lemma E.12b (interpolation errors)
> (i) If $f\in C^2[a,b]$ and $p_1$ is the linear interpolant at $a,b$, then for each $x\in(a,b)$ there is $\eta_x\in(a,b)$ with $f(x)-p_1(x)=\frac12f''(\eta_x)(x-a)(x-b)$.
>
> (ii) If $f\in C^4[a,b]$, there is a unique cubic $p_3$ with $p_3=f$ at $a,m,b$ and $p_3'(m)=f'(m)$. For each $x\in(a,b)\setminus\{m\}$ there is $\eta_x\in(a,b)$ with $f(x)-p_3(x)=\frac1{24}f^{(4)}(\eta_x)(x-a)(x-m)^2(x-b)$.

*Proof.* (ii) *Existence.* Let $p_2$ interpolate $f$ at $a,m,b$, and put $p_3=p_2+c\,(x-a)(x-m)(x-b)$. The cubic term vanishes at the three nodes, and its derivative at $m$ is $(m-a)(m-b)=-h^2\ne0$. So exactly one $c$ makes $p_3'(m)=f'(m)$. Uniqueness follows because the difference of two such cubics would have a double zero at $m$ and simple zeros at $a,b$, which is four zeros counted with multiplicity, so it is identically zero.

*Error.* Fix $x\notin\{a,m,b\}$ and let $\omega(t)=(t-a)(t-m)^2(t-b)$, so $\omega(x)\ne0$. Choose $\kappa$ with $F(t)=f(t)-p_3(t)-\kappa\,\omega(t)$ vanishing at $t=x$. Then $F$ vanishes at the four distinct points $a,m,b,x$, and also $F'(m)=0$. Rolle's theorem on the three gaps between consecutive points gives three zeros of $F'$, each strictly inside a gap and therefore different from $m$. Together with $m$ itself, $F'$ has at least four distinct zeros. Repeating Rolle, $F''$ has three, $F'''$ two, and $F^{(4)}$ at least one, $\eta_x\in(a,b)$. Since $p_3^{(4)}=0$ and $\omega^{(4)}=24$, we get $0=f^{(4)}(\eta_x)-24\kappa$. (i) is the same argument with $\omega(t)=(t-a)(t-b)$, three zeros $a,x,b$, and two differentiations. ∎

**Proof of the theorem**

*Part 1.* By Lemma E.12b(i), $I-Q=\int_a^b(f-p_1)\,dx$, because the trapezoidal rule is the exact integral of $p_1$. Apply Lemma E.12a with $\phi=\frac12f''$, $w=(x-a)(x-b)\le0$ on $[a,b]$, and $q=(f-p_1)/w$ on $(a,b)$. This $q$ is continuous on $(a,b)$, takes values $\frac12f''(\eta_x)$ between $\min\phi$ and $\max\phi$, and $qw=f-p_1$ is integrable. Since $\int_a^b(x-a)(x-b)\,dx=-\frac{(b-a)^3}6$,

$$
I-Q=\tfrac12f''(\xi)\cdot\Big(-\frac{(b-a)^3}6\Big)=-\frac{(b-a)^3}{12}f''(\xi).
$$

If $\xi$ lands on an endpoint, the intermediate value theorem applied on $(a,b)$ to the continuous $f''$ moves it inside.

*Part 2.* **Simpson's rule is exact for cubics.** By linearity it suffices to check the basis $1,\ (x-m),\ (x-m)^2,\ (x-m)^3$. The exact integrals over $[m-h,m+h]$ are $2h,\ 0,\ \frac23h^3,\ 0$. The rule gives $\frac h3(1+4+1)=2h$, $\frac h3(-h+0+h)=0$, $\frac h3(h^2+0+h^2)=\frac23h^3$ and $\frac h3(-h^3+0+h^3)=0$. **Reduction.** The rule uses only $f(a),f(m),f(b)$, where $f=p_3$. So $Q(f)=Q(p_3)=\int_a^bp_3\,dx$, and $I-Q=\int_a^b(f-p_3)\,dx$. **Weighted mean value.** Lemma E.12a applies with $\phi=f^{(4)}/24$ and $w=(x-a)(x-m)^2(x-b)\le0$ on $[a,b]$. With $x=m+t$,

$$
\int_a^bw\,dx=\int_{-h}^{h}(t+h)\,t^2\,(t-h)\,dt=\int_{-h}^{h}\big(t^4-h^2t^2\big)\,dt=\frac{2h^5}5-\frac{2h^5}3=-\frac{4h^5}{15},
$$

so $I-Q=\frac{f^{(4)}(\xi)}{24}\cdot\big(-\frac{4h^5}{15}\big)=-\frac{h^5}{90}f^{(4)}(\xi)$.

*Part 3.* The composite trapezoidal error is the sum of the $n$ panel errors, each $-\frac{h^3}{12}f''(\xi_i)$ by part 1. The average $\frac1n\sum_if''(\xi_i)$ lies between $\min f''$ and $\max f''$, so it equals $f''(\xi)$ for some $\xi$ by the intermediate value theorem. The total is $-\frac{nh^3}{12}f''(\xi)=-\frac{(b-a)h^2}{12}f''(\xi)$, since $nh=b-a$. For Simpson, the $n/2$ double panels each contribute $-\frac{h^5}{90}f^{(4)}(\xi_i)$ by part 2, and the same averaging gives $-\frac n2\cdot\frac{h^5}{90}f^{(4)}(\xi')=-\frac{(b-a)h^4}{180}f^{(4)}(\xi')$. ∎

```desmos-graph
left=-1.2; right=1.2; bottom=-1.1; top=0.4
---
y=(x+1)(x-1)|blue|-1<x<1
y=x^{2}(x^{2}-1)|red|-1<x<1
y=0|black
```

*The error kernels on $[-1,1]$: $(x-a)(x-b)$ for the trapezoidal rule (blue) and $(x-a)(x-m)^2(x-b)$ for Simpson's rule (red). Neither changes sign, which is exactly what Lemma E.12a needs. Their integrals, $-\frac43$ and $-\frac4{15}$, produce the constants $\frac1{12}$ and $\frac1{90}$.*

**Remarks**

- *Numerical check.* For $\int_0^1e^x\,dx=e-1$, the single trapezoid has error $-0.14086=-\frac1{12}e^\xi$, giving $\xi=0.525\in(0,1)$. Single-panel Simpson has error $-5.79\times10^{-4}=-\frac{(0.5)^5}{90}e^\xi$, giving $\xi=0.512$. Both are consistent with the theorem.
- *Counterexample: smoothness matters.* For $f(x)=\sqrt x$ on $[0,1]$, $f''$ is unbounded and part 3 does not apply. The composite trapezoidal errors for $n=8,16,32,64$ are $8.5\times10^{-3},\ 3.1\times10^{-3},\ 1.1\times10^{-3},\ 4.0\times10^{-4}$. They shrink by about $2.8\approx2^{1.5}$ per halving, an $O(h^{1.5})$ rate instead of $O(h^2)$. A singularity at an endpoint degrades the order. Subtracting the singularity, changing variables, or Gaussian quadrature restores it ([Section 22.8](../part-7-numerical-computational/22-numerical-methods.md)).
- *The bound is not always sharp.* For smooth periodic integrands over a full period, the trapezoidal rule converges far faster than $h^2$. Every term of the Euler–Maclaurin expansion vanishes ([Section 22.8](../part-7-numerical-computational/22-numerical-methods.md)). The theorem gives a guarantee, not the true error in every case.
- *Where used.* [Section 22.8](../part-7-numerical-computational/22-numerical-methods.md) (Newton–Cotes and composite rules, Richardson/Romberg). The same Rolle-counting proof gives the general interpolation error formula of [Section 22.5](../part-7-numerical-computational/22-numerical-methods.md).

## E.13 Euler's method: global $O(h)$ convergence and absolute stability

Consider the initial-value problem $y'=f(t,y)$, $y(t_0)=y_0$, on $[t_0,T]$. Euler's method with step $h>0$ is $y_{n+1}=y_n+hf(t_n,y_n)$ with $t_n=t_0+nh$ ([Section 22.9](../part-7-numerical-computational/22-numerical-methods.md)). The global error is $e_n=y(t_n)-y_n$.

**Statement**

> [!theorem] Convergence of Euler's method
> Assume that $f$ is continuous on $[t_0,T]\times\mathbb R$ and **Lipschitz** in $y$, i.e. $\lvert f(t,u)-f(t,v)\rvert\le L\lvert u-v\rvert$ for all $t,u,v$, with $L>0$. Assume also that the exact solution $y$ is twice continuously differentiable with $\lvert y''\rvert\le M$ on $[t_0,T]$. If the method is started from $y_0$ with initial error $e_0$, then for every $n$ with $t_n\le T$
> $$
> \lvert e_n\rvert\le e^{L(t_n-t_0)}\lvert e_0\rvert+\frac{Mh}{2L}\Big(e^{L(t_n-t_0)}-1\Big).
> $$
> In particular, with $e_0=0$, $\max_{t_n\le T}\lvert e_n\rvert\le C\,h$ with $C=\frac M{2L}(e^{L(T-t_0)}-1)$ independent of $h$. Euler's method is convergent of order one.

> [!theorem] Lemma E.13a (discrete Gronwall inequality)
> If $a_n\ge0$ and $a_{n+1}\le(1+hL)a_n+B$ for $n\ge0$, with $h,L>0$ and $B\ge0$, then
> $$
> a_n\le(1+hL)^na_0+B\,\frac{(1+hL)^n-1}{hL}\le e^{nhL}a_0+\frac{B}{hL}\big(e^{nhL}-1\big).
> $$

*Proof.* Induction gives $a_n\le(1+hL)^na_0+B\sum_{k=0}^{n-1}(1+hL)^k$. The geometric sum equals $\big((1+hL)^n-1\big)/(hL)$. Finally, $1+x\le e^x$ for all real $x$, because $e^x$ is convex and $1+x$ is its tangent at $0$. Hence $(1+hL)^n\le e^{nhL}$. ∎

**Proof of the theorem**

*Step 1 (local truncation error).* By Taylor's theorem with Lagrange remainder and $y'(t_n)=f(t_n,y(t_n))$,

$$
y(t_{n+1})=y(t_n)+hf\big(t_n,y(t_n)\big)+\tau_{n+1},\qquad\tau_{n+1}=\tfrac{h^2}2y''(\xi_n),\qquad\lvert\tau_{n+1}\rvert\le\tfrac M2h^2 .
$$

*Step 2 (error recursion).* Subtract the Euler step $y_{n+1}=y_n+hf(t_n,y_n)$:

$$
e_{n+1}=e_n+h\big[f(t_n,y(t_n))-f(t_n,y_n)\big]+\tau_{n+1}.
$$

By the Lipschitz condition, the bracket is at most $L\lvert e_n\rvert$ in modulus. So $\lvert e_{n+1}\rvert\le(1+hL)\lvert e_n\rvert+\frac M2h^2$.

*Step 3 (accumulate).* Apply Lemma E.13a with $a_n=\lvert e_n\rvert$ and $B=\frac M2h^2$. Then $B/(hL)=Mh/(2L)$, and $nh=t_n-t_0$. ∎

The structure of the proof is the general principle "consistency + stability ⟹ convergence" of [Section 22.11](../part-7-numerical-computational/22-numerical-methods.md). There are $O(1/h)$ local errors, each $O(h^2)$. Each is amplified by at most $e^{L(T-t_0)}$ by the time it reaches $T$.

**Statement (absolute stability)**

> [!theorem] Absolute stability of Euler's method
> 1. For the test equation $y'=\lambda y$, $\lambda\in\mathbb C$, Euler's method gives $y_n=(1+h\lambda)^ny_0$. The numerical solution is bounded for every $y_0$ iff $\lvert1+h\lambda\rvert\le1$, and tends to $0$ for every $y_0$ iff $\lvert1+h\lambda\rvert<1$. The region of absolute stability $\{z=h\lambda:\lvert1+z\rvert\le1\}$ is the closed disc of radius $1$ centred at $-1$.
> 2. For real $\lambda<0$ the condition is $0<h\le2/\lvert\lambda\rvert$. The numerical solution decays without sign changes iff $h\le1/\lvert\lambda\rvert$.
> 3. For $\mathbf y'=\mathbf A\mathbf y$ with $\mathbf A$ diagonalisable, with eigenvalues $\lambda_1,\dots,\lambda_N$, the Euler iterates $\mathbf y_n=(\mathbf I+h\mathbf A)^n\mathbf y_0$ are bounded for every $\mathbf y_0$ iff $\lvert1+h\lambda_i\rvert\le1$ for every $i$.

**Proof**

*Part 1.* With $f(t,y)=\lambda y$ the Euler step is $y_{n+1}=(1+h\lambda)y_n$, so by induction $y_n=(1+h\lambda)^ny_0$ and $\lvert y_n\rvert=\lvert1+h\lambda\rvert^n\lvert y_0\rvert$. For $y_0\ne0$, a geometric sequence $\rho^n$ with $\rho\ge0$ is bounded iff $\rho\le1$ and tends to $0$ iff $\rho<1$. Writing $z=h\lambda=x+jy$, the condition $\lvert1+z\rvert\le1$ reads $(x+1)^2+y^2\le1$.

*Part 2.* For real $z=h\lambda<0$, $\lvert1+z\rvert\le1$ means $-1\le1+z\le1$, i.e. $-2\le z\le0$, i.e. $h\le2/\lvert\lambda\rvert$. The factor $1+z$ is non-negative, so the iterates keep the sign of $y_0$, iff $z\ge-1$.

*Part 3.* Write $\mathbf A=\mathbf V\boldsymbol\Lambda\mathbf V^{-1}$ with $\mathbf V$ invertible. Then $\mathbf I+h\mathbf A=\mathbf V(\mathbf I+h\boldsymbol\Lambda)\mathbf V^{-1}$ and $(\mathbf I+h\mathbf A)^n=\mathbf V(\mathbf I+h\boldsymbol\Lambda)^n\mathbf V^{-1}$. In the coordinates $\mathbf w=\mathbf V^{-1}\mathbf y$, the iteration decouples into $w_{i,n}=(1+h\lambda_i)^nw_{i,0}$, and part 1 applies to each component. Since $\lVert\mathbf V^{-1}\rVert^{-1}\lVert\mathbf y_n\rVert\le\lVert\mathbf w_n\rVert\le\lVert\mathbf V^{-1}\rVert\lVert\mathbf y_n\rVert$, the sequence $\mathbf y_n$ is bounded iff $\mathbf w_n$ is. ∎

```desmos-graph
left=-3.2; right=1; bottom=-1.6; top=1.6
---
(x+1)^{2}+y^{2}\le1|blue
x=0|black|dashed
(-0.5,0)|blue|label:z=-0.5
(-1.5,0)|orange|label:z=-1.5
(-2.5,0)|red|label:z=-2.5
```

*Region of absolute stability of explicit Euler in the complex $z=h\lambda$ plane: the disc $\lvert1+z\rvert\le1$. The three points are the step sizes of the chapter's figure for $y'=-2y$. $z=-0.5$ gives monotone decay, $z=-1.5$ oscillating decay, and $z=-2.5$ lies outside the disc, so the iterates grow.*

**Remarks**

- *Numerical check of the convergence bound.* For $y'=y$, $y(0)=1$ on $[0,1]$, take $L=1$ and $M=e$. At $t=1$ the bound is $\frac{eh}2(e-1)=2.335h$. The actual errors $e-(1+h)^{1/h}$ are $0.1245,\ 0.0650,\ 0.01347$ for $h=0.1,\ 0.05,\ 0.01$. Each is below the bound, and the ratio error$/h$ tends to $e/2=1.359$, confirming first order.
- *Counterexample: the Lipschitz hypothesis.* $y'=\sqrt{\lvert y\rvert}$, $y(0)=0$, is not Lipschitz at $y=0$. It has the solutions $y\equiv0$ and $y=t^2/4$ (among others). Euler from $y_0=0$ produces $y_n\equiv0$ and converges to one solution. If $y_0=\varepsilon>0$, it follows the other branch instead. Without uniqueness, "the error" is not even well defined.
- *Convergence versus stability.* The convergence theorem is about $h\to0$ at fixed $T$. Its constant $e^{L(T-t_0)}$ can be astronomically large. The stability theorem is about fixed $h$ and $n\to\infty$. For a stiff system with $\lambda=-1000$, part 2 forces $h\le2\times10^{-3}$ for as long as the integration runs, even after the fast mode has decayed. That is the motivation for implicit methods ([Section 22.9](../part-7-numerical-computational/22-numerical-methods.md), "Stiffness").
- *Counterexample: diagonalisability in part 3.* For $\mathbf A=\begin{pmatrix}\lambda&1\\0&\lambda\end{pmatrix}$ with $\lvert1+h\lambda\rvert=1$, $(\mathbf I+h\mathbf A)^n=\begin{pmatrix}r^n&nhr^{n-1}\\0&r^n\end{pmatrix}$ with $r=1+h\lambda$. It grows linearly in $n$ although every eigenvalue lies on the boundary of the stability disc. Strict inequality $\lvert1+h\lambda_i\rvert<1$ suffices in general.
- *Rounding.* If each step also commits a rounding error of size $\epsilon$, then $B=\frac M2h^2+\epsilon$ in Step 3, and the bound becomes $\big(\frac{Mh}2+\frac\epsilon h\big)\frac{e^{L(T-t_0)}-1}L$. It is minimised near $h=\sqrt{2\epsilon/M}$. Making $h$ smaller beyond this point makes the answer worse.
- *Where used.* Euler's method and absolute stability in [Section 22.9](../part-7-numerical-computational/22-numerical-methods.md), zero-stability and the Lax equivalence theorem in [Section 22.11](../part-7-numerical-computational/22-numerical-methods.md), and solver choice in [Section 25.9](../part-7-numerical-computational/25-mathematical-modeling.md).

## E.14 The Monte Carlo $1/\sqrt N$ error rate

**Statement**

> [!theorem] Error of the Monte Carlo estimator
> Let $\mathbf X_1,\dots,\mathbf X_N$ be i.i.d. copies of a random vector $\mathbf X$, and let $g$ be a function with $\mu=E[g(\mathbf X)]$ and $\sigma^2=\operatorname{Var}g(\mathbf X)<\infty$. Let
> $$
> \hat\mu_N=\frac1N\sum_{i=1}^Ng(\mathbf X_i),\qquad s_N^2=\frac1{N-1}\sum_{i=1}^N\big(g(\mathbf X_i)-\hat\mu_N\big)^2 .
> $$
> Then:
> 1. $E[\hat\mu_N]=\mu$, and the root-mean-square error is **exactly** $\sqrt{E[(\hat\mu_N-\mu)^2]}=\sigma/\sqrt N$;
> 2. for every $N$ and every $k>0$, $P\big(\lvert\hat\mu_N-\mu\rvert\ge k\sigma/\sqrt N\big)\le1/k^2$;
> 3. if $\sigma>0$, then for every $z>0$, $P\big(\lvert\hat\mu_N-\mu\rvert\le z\sigma/\sqrt N\big)\to2\Phi(z)-1$ as $N\to\infty$;
> 4. $E[s_N^2]=\sigma^2$ and $s_N\to\sigma$ in probability. Consequently, part 3 also holds with $\sigma$ replaced by $s_N$, and $\hat\mu_N\pm1.96\,s_N/\sqrt N$ is an asymptotic $95\%$ confidence interval.

**Proof**

*Part 1.* The $Y_i=g(\mathbf X_i)$ are i.i.d. with mean $\mu$ and variance $\sigma^2$. By E.4 (part 1 of the proof), $E[\hat\mu_N]=\mu$ and $\operatorname{Var}\hat\mu_N=\sigma^2/N$. Independence implies uncorrelatedness, so the variance identity is exact, not just a bound. Since $E[\hat\mu_N]=\mu$, the mean-square error equals the variance.

*Part 2.* Chebyshev's inequality (E.3) applied to $\hat\mu_N$, whose standard deviation is $\sigma/\sqrt N$.

*Part 3.* By the central limit theorem (E.5) applied to the $Y_i$, $Z_N=\sqrt N(\hat\mu_N-\mu)/\sigma\to\mathcal N(0,1)$ in distribution. Because $\Phi$ is continuous, $P(-z\le Z_N\le z)=F_{Z_N}(z)-F_{Z_N}(-z)+P(Z_N=-z)$. The first two terms converge to $\Phi(z)-\Phi(-z)$. For the atom term, $0\le P(Z_N=-z)\le F_{Z_N}(-z)-F_{Z_N}(-z-\eta)\to\Phi(-z)-\Phi(-z-\eta)$ for every $\eta>0$, which can be made arbitrarily small. So $P(\lvert Z_N\rvert\le z)\to2\Phi(z)-1$.

*Part 4.* Unbiasedness is E.6. For consistency, expand as in E.6, Step 1:

$$
s_N^2=\frac N{N-1}\Big[\frac1N\sum_{i=1}^NY_i^2-\hat\mu_N^2\Big].
$$

By Khinchin's weak law (E.4, part 2), which applies because $E\lvert Y_i^2\rvert=\sigma^2+\mu^2<\infty$, we have $\frac1N\sum Y_i^2\to\sigma^2+\mu^2$ and $\hat\mu_N\to\mu$ in probability. Continuous functions preserve convergence in probability (the continuous mapping theorem, quoted). So $s_N^2\to(\sigma^2+\mu^2)-\mu^2=\sigma^2$, and $s_N=\sqrt{s_N^2}\to\sigma$. Slutsky's theorem (quoted: if $Z_N\to Z$ in distribution and $C_N\to c\ne0$ in probability, then $Z_N/C_N\to Z/c$ in distribution) applied with $C_N=s_N/\sigma\to1$ shows that $\sqrt N(\hat\mu_N-\mu)/s_N\to\mathcal N(0,1)$. The interval statement follows as in part 3. ∎

**Remarks**

- *The rate is exact.* Part 1 is an identity, not a bound. For plain Monte Carlo with i.i.d. samples, the RMS error is $\sigma N^{-1/2}$ for every $N$ and cannot be improved by any choice of $N$. Halving the error costs four times the samples, and one more decimal digit costs a hundred times. The only levers are $\sigma$ (variance reduction, which replaces $g$ by an estimator with the same mean and smaller variance) and abandoning independence (quasi-Monte Carlo) ([Section 24.6](../part-7-numerical-computational/24-computational-mathematics.md)).
- *No dimension in the rate.* Neither the dimension of $\mathbf X$ nor the smoothness of $g$ appears in the theorem. Only $\sigma$ does. By contrast, a product trapezoidal rule in $d$ dimensions has error $O(N^{-2/d})$ (E.12, part 3, in each coordinate). This is why Monte Carlo wins for $d>4$.
- *Numerical check.* For $\int_0^1e^x\,dx$ with uniform samples, $\sigma=\sqrt{\frac12(e^2-1)-(e-1)^2}=0.49197$. Over $2000$ independent repetitions, the observed RMS errors were $0.0496,\ 0.0152,\ 0.00481$ for $N=10^2,\ 10^3,\ 10^4$. The predicted values $\sigma/\sqrt N$ are $0.0492,\ 0.0156,\ 0.00492$.
- *Counterexample: infinite variance.* For $g(x)=x^{-1/2}$ with $X\sim\mathcal U(0,1)$, $\mu=2$ is finite but $E[g^2]=\int_0^1x^{-1}dx=\infty$. Khinchin's law still gives $\hat\mu_N\to2$ in probability, but parts 1–4 fail. The error decays more slowly than $N^{-1/2}$ (here like $N^{-1/2}\sqrt{\log N}$), and $s_N$ systematically misreports it. Importance sampling with a density $\propto x^{-1/2}$ removes the singularity and gives zero variance.
- *Chebyshev versus CLT.* Part 2 holds for every $N$ but is conservative. With $k=4.47$ it guarantees $95\%$ coverage. Part 3 needs only $z=1.96$, but it is asymptotic, and for skewed integrands with small $N$ the actual coverage can fall short (Berry–Esseen, E.5).
- *Where used.* [Section 24.6](../part-7-numerical-computational/24-computational-mathematics.md) (the basic estimator and its error bar), Monte Carlo tolerance analysis and uncertainty propagation ([Sections 24.6 and 24.8](../part-7-numerical-computational/24-computational-mathematics.md)), and GUM Supplement 1 in [Section 20.10](../part-6-probability-statistics/20-statistical-analysis.md).

## Further reading

- P. Billingsley, *Probability and Measure*, 3rd ed., Wiley, 1995. Measure-theoretic foundations, Lévy's continuity theorem, and the Lindeberg CLT.
- R. Durrett, *Probability: Theory and Examples*, 5th ed., Cambridge University Press, 2019. Laws of large numbers, characteristic functions, and the CLT with full proofs.
- G. Grimmett and D. Stirzaker, *Probability and Random Processes*, 4th ed., Oxford University Press, 2020. An accessible rigorous treatment, including stationary processes and spectral theory.
- A. Papoulis and S. U. Pillai, *Probability, Random Variables and Stochastic Processes*, 4th ed., McGraw-Hill, 2002. The engineering standard for Wiener–Khinchin and linear systems with random inputs.
- G. A. F. Seber and A. J. Lee, *Linear Regression Analysis*, 2nd ed., Wiley, 2003. Least squares, Gauss–Markov and the distribution theory of regression.
- JCGM 100:2008, *Evaluation of Measurement Data — Guide to the Expression of Uncertainty in Measurement* (GUM), and JCGM 101:2008 (Supplement 1, Monte Carlo propagation).
- E. Süli and D. Mayers, *An Introduction to Numerical Analysis*, Cambridge University Press, 2003. Rigorous proofs for root finding, interpolation, quadrature and ODE methods.
- E. Hairer, S. P. Nørsett and G. Wanner, *Solving Ordinary Differential Equations I*, 2nd ed., Springer, 1993. Convergence and stability theory of one-step methods.
- A. B. Owen, *Monte Carlo Theory, Methods and Examples*, 2013 (online). Monte Carlo error analysis, variance reduction and quasi-Monte Carlo.
