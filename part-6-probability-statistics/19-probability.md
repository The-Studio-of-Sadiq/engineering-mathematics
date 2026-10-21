# Chapter 19 — Probability

[← Previous](../part-5-transforms/18-complex-analysis.md) · [Table of Contents](../00-front/table-of-contents.md) · [Next →](20-statistical-analysis.md)

> **Prerequisites:** Ch. 1 (sets, functions, sums), Ch. 5 (integral calculus), Ch. 6 (multiple integrals), Ch. 16 (Fourier transforms, for the characteristic function)
> **Learning outcomes:** After this chapter you can …
> - build a probability model — sample space, events, probability measure — from the axioms, and count outcomes with permutations and combinations;
> - compute conditional probabilities, apply the law of total probability, and invert conditioning with Bayes' theorem (binary channels, fault detectors, MAP decisions);
> - decide whether events are independent and compute the reliability of series, parallel, $k$-out-of-$n$ and bridge systems;
> - describe discrete, continuous and mixed random variables by their PMF, PDF and CDF, and derive the distribution of a function of a random variable;
> - compute expectations, variances and higher moments, use the Markov and Chebyshev bounds, and obtain moments from the moment-generating and characteristic functions;
> - set up a joint distribution of two random variables, find marginals and conditionals, and compute the distribution of a sum of independent random variables.

Every measurement an engineer makes is uncertain. The voltage across a resistor fluctuates because electrons are in thermal motion; a bit sent over a wireless link is sometimes received in error; a batch of capacitors has values scattered around the nominal; a power-supply component fails after an unpredictable number of hours; packets arrive at a router at irregular instants. None of these can be predicted exactly, but all of them are *regular in the aggregate*: the fraction of bits in error over a long run settles to a stable number, the histogram of capacitor values has a repeatable shape, the average lifetime of a large population of components is well defined. **Probability theory** is the mathematics of that aggregate regularity. It lets us replace "we cannot know" by quantitative statements — "the bit-error rate is $10^{-6}$", "the system survives five years with probability $0.98$", "the noise exceeds the threshold $0.3\,\%$ of the time" — on which designs can be based.

This chapter builds the theory from the ground up. Sections 19.1–19.4 work with *events*: the axioms, counting, conditioning, Bayes' theorem and independence. Sections 19.5–19.8 introduce *random variables* — numerical outcomes — and the three functions that describe them (PMF, PDF, CDF), together with the most important families of distributions and the first ideas of joint distributions. Section 19.9 condenses a distribution into a few numbers: the mean, the variance, the higher moments, and the transforms (moment-generating and characteristic functions) that generate them. [Ch. 20](20-statistical-analysis.md) turns the theory around — from a model to data, and from data back to a model — and catalogues the common distributions in more depth; [Ch. 21](21-random-processes.md) lets the random variable depend on time.

Throughout, outcomes are denoted $\zeta$, events by capital letters $A, B, \dots$, random variables by capital letters $X, Y, \dots$ and their values by the corresponding lowercase letters $x, y, \dots$.

## 19.1 Sample spaces and events

### Motivation

Before we can attach numbers to uncertainty we must say precisely *what* is uncertain. "The probability that the receiver makes an error" is meaningless until we specify the experiment (one bit? one packet? one hour of operation?), list what can happen, and agree which collections of outcomes we want to talk about. This bookkeeping is the sample space and its events. Getting it right removes most of the paradoxes and wrong answers that plague informal probability arguments.

### Experiments, outcomes and sample spaces

> [!definition] Random experiment, outcome, sample space
> A **random experiment** is a procedure whose result cannot be predicted with certainty but whose set of possible results is known in advance. Each possible result is an **outcome** $\zeta$. The set of all outcomes is the **sample space** $\Omega$ (some texts write $S$).

The outcomes must be **mutually exclusive** (exactly one occurs on each trial) and **collectively exhaustive** (one of them always occurs). How finely to describe the outcomes is a modelling choice: for a single coin toss $\Omega=\{H,T\}$ suffices, but if you care about where the coin lands you need a richer description.

Sample spaces come in three kinds:

| Kind | Example experiment | Sample space |
|---|---|---|
| Finite | Send a 3-bit word and record which bits are in error | $\Omega=\{000,001,\dots,111\}$, 8 outcomes (1 = bit in error) |
| Countably infinite | Retransmit a packet until it is received correctly; record the number of attempts | $\Omega=\{1,2,3,\dots\}$ |
| Uncountable (continuum) | Measure the lifetime of a capacitor in hours | $\Omega=[0,\infty)$ |
| Uncountable, multidimensional | Sample the in-phase and quadrature noise voltages of a receiver | $\Omega=\mathbb{R}^2$ |

### Events and set operations

> [!definition] Event
> An **event** is a subset $A\subseteq\Omega$ of outcomes. We say the event $A$ **occurs** on a trial if the outcome $\zeta$ of that trial satisfies $\zeta\in A$. The whole space $\Omega$ is the **certain event**; the empty set $\varnothing$ is the **impossible event**; a single outcome $\{\zeta\}$ is an **elementary event**.

Because events are sets, the algebra of sets ([Ch. 1](../part-1-foundations/01-mathematical-language-and-algebra.md)) becomes the logic of events:

| Set notation | Read as | Event occurs when … |
|---|---|---|
| $A\cup B$ | $A$ or $B$ (union) | at least one of $A$, $B$ occurs |
| $A\cap B$ (also $AB$) | $A$ and $B$ (intersection) | both occur |
| $A^c$ (also $\bar A$) | not $A$ (complement) | $A$ does not occur |
| $A\setminus B = A\cap B^c$ | $A$ but not $B$ | $A$ occurs and $B$ does not |
| $A\subseteq B$ | $A$ implies $B$ | whenever $A$ occurs, $B$ occurs |
| $A\cap B=\varnothing$ | $A$, $B$ mutually exclusive (disjoint) | they never occur together |

The rules that matter most in calculation are the **distributive laws** and **De Morgan's laws**:

$$
A\cap(B\cup C)=(A\cap B)\cup(A\cap C),\qquad (A\cup B)^c=A^c\cap B^c,\qquad (A\cap B)^c=A^c\cup B^c .
$$

De Morgan's laws extend to any number of events: "not (at least one occurs)" is "none occurs", $\left(\bigcup_i A_i\right)^c=\bigcap_i A_i^c$. This is the single most useful identity in reliability and error analysis, because "at least one component fails" is usually hard to compute directly while "no component fails" is easy.

A collection of events $B_1,B_2,\dots,B_n$ is a **partition** of $\Omega$ if the events are pairwise disjoint and their union is $\Omega$. Exactly one event of a partition occurs on each trial. Partitions are the scaffolding of the law of total probability (Section 19.2).

> [!example] Worked Example 19.1 — Events in a 3-bit word
> A 3-bit word is transmitted; outcome $\zeta=b_1b_2b_3$ with $b_i=1$ if bit $i$ is received in error. Let $A$ = "bit 1 is in error", $B$ = "exactly one bit is in error", $C$ = "the word is received correctly".
>
> **Step 1 — list the sets.** $A=\{100,101,110,111\}$, $B=\{100,010,001\}$, $C=\{000\}$.
>
> **Step 2 — combine.** $A\cap B=\{100\}$ ("only bit 1 is wrong"); $A\cup B=\{100,101,110,111,010,001\}$; $C^c$ = "at least one error" $=\Omega\setminus\{000\}$, which contains 7 outcomes.
>
> **Step 3 — check disjointness.** $B\cap C=\varnothing$: a word cannot have exactly one error and no error. $A\cap C=\varnothing$ as well. The events "0 errors", "1 error", "2 errors", "3 errors" form a partition of $\Omega$ into sets of sizes $1,3,3,1$.

### Which subsets are events? ($\sigma$-algebras)

For finite and countable $\Omega$ every subset can be an event. For a continuum such as $\Omega=[0,1]$ it turns out to be impossible to assign a probability to *every* subset consistently with the axioms below (a classical result of measure theory). We therefore restrict attention to a family $\mathcal{F}$ of subsets, the **events**, that is closed under the operations we need.

> [!definition] $\sigma$-algebra (event field)
> A family $\mathcal{F}$ of subsets of $\Omega$ is a **$\sigma$-algebra** if
> 1. $\Omega\in\mathcal{F}$;
> 2. $A\in\mathcal{F}\Rightarrow A^c\in\mathcal{F}$;
> 3. $A_1,A_2,\dots\in\mathcal{F}\Rightarrow \bigcup_{i=1}^{\infty}A_i\in\mathcal{F}$.
>
> By De Morgan's laws it is then also closed under countable intersections. On $\mathbb{R}$ the standard choice is the **Borel** $\sigma$-algebra, the smallest one containing all intervals; it contains every set an engineer will ever write down.

For practical work you may treat "event" as "any subset you can describe". The $\sigma$-algebra is there to guarantee that limits of events — "the error count *eventually* stays below 5", "the voltage *ever* exceeds 1 V" — are themselves events.

### The axioms of probability

Probability was put on a rigorous footing by Kolmogorov (1933), who took the following three properties as axioms and derived everything else from them.

> [!definition] Probability space (Kolmogorov axioms)
> A **probability space** is a triple $(\Omega,\mathcal{F},P)$ where $P$ assigns to every event $A\in\mathcal{F}$ a real number $P(A)$ such that
> 1. **Non-negativity:** $P(A)\ge 0$;
> 2. **Normalization:** $P(\Omega)=1$;
> 3. **Countable additivity:** if $A_1,A_2,\dots$ are pairwise disjoint, then
>
> $$
> P\!\left(\bigcup_{i=1}^{\infty}A_i\right)=\sum_{i=1}^{\infty}P(A_i).
> $$

The axioms do not say *what* number to assign to an event; they only constrain the assignment to be internally consistent. The numbers come from modelling: symmetry ("all 36 outcomes of two dice are equally likely"), physics (Boltzmann statistics of thermal noise), or measurement (a long-run relative frequency of failures). The **relative-frequency interpretation** — if the experiment is repeated $N$ times and $A$ occurs $N_A$ times, then $N_A/N\to P(A)$ as $N\to\infty$ — is what connects the axioms to the laboratory, and it is a *theorem* of the theory (the law of large numbers, Section 19.9 and [Ch. 20](20-statistical-analysis.md)) rather than a definition.

### Consequences of the axioms

Every rule of probability you will use follows from the three axioms by short arguments.

> [!theorem] Elementary properties
> For events $A,B$ in a probability space:
> 1. $P(\varnothing)=0$.
> 2. $P(A^c)=1-P(A)$.
> 3. If $A\subseteq B$ then $P(A)\le P(B)$ (monotonicity); in particular $0\le P(A)\le 1$.
> 4. $P(A\cup B)=P(A)+P(B)-P(A\cap B)$ (addition rule).
> 5. $P(A\cup B)\le P(A)+P(B)$, and more generally $P\left(\bigcup_i A_i\right)\le\sum_i P(A_i)$ (union, or Boole, bound).

*Proof sketch.* (1) Take $A_1=\Omega$ and $A_i=\varnothing$ for $i\ge2$; additivity gives $P(\Omega)=P(\Omega)+\sum_{i\ge2}P(\varnothing)$, forcing $P(\varnothing)=0$. (2) $A$ and $A^c$ are disjoint with union $\Omega$, so $P(A)+P(A^c)=1$. (3) Write $B=A\cup(B\setminus A)$, a disjoint union, so $P(B)=P(A)+P(B\setminus A)\ge P(A)$. (4) Split $A\cup B$ into the three disjoint pieces $A\setminus B$, $A\cap B$, $B\setminus A$; also $P(A)=P(A\setminus B)+P(A\cap B)$ and similarly for $B$. Adding the last two and subtracting gives the result: the intersection has been counted twice. (5) follows from (4) because $P(A\cap B)\ge 0$, and by induction for finitely many events. $\square$

The addition rule generalizes to the **inclusion–exclusion principle**:

$$
P\!\left(\bigcup_{i=1}^{n}A_i\right)=\sum_i P(A_i)-\sum_{i<j}P(A_iA_j)+\sum_{i<j<k}P(A_iA_jA_k)-\cdots+(-1)^{n+1}P(A_1A_2\cdots A_n).
$$

For three events, $P(A\cup B\cup C)=P(A)+P(B)+P(C)-P(AB)-P(AC)-P(BC)+P(ABC)$. Truncating the series after the first sum gives the union bound (an upper bound); truncating after the second sum gives a lower bound (the **Bonferroni inequalities** alternate in this way).

A final consequence is **continuity of probability**: if $A_1\subseteq A_2\subseteq\cdots$ then $P\left(\bigcup_n A_n\right)=\lim_n P(A_n)$, and if $A_1\supseteq A_2\supseteq\cdots$ then $P\left(\bigcap_n A_n\right)=\lim_n P(A_n)$. This is what will make the CDF right-continuous in Section 19.8.

> [!example] Worked Example 19.2 — Two test stations (addition rule)
> Boards from a production line pass through an in-circuit test (ICT) and a functional test (FT). Records show $P(\text{fails ICT})=0.05$, $P(\text{fails FT})=0.04$ and $P(\text{fails both})=0.01$. Find the probability that a board fails at least one test, fails exactly one test, and passes both.
>
> **Step 1 — at least one.** Let $I$ and $F$ be the failure events. By the addition rule,
> $$
> P(I\cup F)=0.05+0.04-0.01=0.08 .
> $$
>
> **Step 2 — exactly one.** "Exactly one" is $(I\cup F)\setminus(I\cap F)$, so its probability is $0.08-0.01=0.07$. Check: $P(I\setminus F)=0.05-0.01=0.04$ and $P(F\setminus I)=0.04-0.01=0.03$; total $0.07$. ✓
>
> **Step 3 — passes both.** By De Morgan, "passes both" is $I^c\cap F^c=(I\cup F)^c$, probability $1-0.08=0.92$.
>
> The union bound would have given $P(I\cup F)\le 0.09$: close, because the overlap is small. When failures are rare, the union bound is nearly exact — the basis of many quick error-rate estimates.

### Equally likely outcomes and counting

When $\Omega$ is finite with $N$ outcomes and symmetry makes them equally likely, additivity forces each outcome to have probability $1/N$, and

$$
P(A)=\frac{\lvert A\rvert}{\lvert\Omega\rvert}=\frac{\text{number of outcomes in }A}{\text{total number of outcomes}} .
$$

Probability then reduces to **counting**, and counting large sets needs a few systematic tools — the subject of combinatorics.

> [!theorem] Multiplication (product) principle
> If a procedure consists of $k$ stages, and stage $i$ can be carried out in $n_i$ ways regardless of how the earlier stages were done, then the whole procedure can be carried out in $n_1n_2\cdots n_k$ ways.

The four basic counts follow by choosing $k$ items from $n$ distinct items, distinguishing whether order matters and whether repetition is allowed:

| | Order matters | Order does not matter |
|---|---|---|
| **With replacement** | $n^k$ | $\dbinom{n+k-1}{k}$ |
| **Without replacement** | $P(n,k)=\dfrac{n!}{(n-k)!}$ | $\dbinom{n}{k}=\dfrac{n!}{k!\,(n-k)!}$ |

*Derivations.* Ordered with replacement: $k$ stages, $n$ choices each, giving $n^k$ (the number of $k$-symbol words over an $n$-letter alphabet; $2^k$ binary words of length $k$). Ordered without replacement: $n$ choices for the first item, $n-1$ for the second, …, $n-k+1$ for the $k$-th, product $n!/(n-k)!$; with $k=n$ this is $n!$, the number of **permutations** of $n$ items. Unordered without replacement: each unordered set of $k$ items corresponds to exactly $k!$ ordered lists, so the number of **combinations** is $P(n,k)/k!=\binom{n}{k}$. Unordered with replacement ("stars and bars"): a multiset of size $k$ from $n$ types is a row of $k$ stars separated into $n$ bins by $n-1$ bars; choosing the positions of the $k$ stars among $n+k-1$ symbols gives $\binom{n+k-1}{k}$.

The binomial coefficient obeys $\binom{n}{k}=\binom{n}{n-k}$, Pascal's rule $\binom{n}{k}=\binom{n-1}{k-1}+\binom{n-1}{k}$, and the **binomial theorem** $(a+b)^n=\sum_{k=0}^{n}\binom{n}{k}a^kb^{n-k}$, which will normalize the binomial distribution in Section 19.6.

> [!definition] Multinomial coefficient
> The number of ways to split $n$ distinct items into groups of sizes $k_1,\dots,k_r$ with $k_1+\cdots+k_r=n$ (equivalently, the number of distinct arrangements of a word with $k_i$ copies of letter $i$) is
> $$
> \binom{n}{k_1,k_2,\dots,k_r}=\frac{n!}{k_1!\,k_2!\cdots k_r!}.
> $$

For example, the number of distinct orders in which 4 identical resistors, 3 identical capacitors and 2 identical inductors can be placed in 9 slots is $9!/(4!\,3!\,2!)=1260$.

For large $n$, factorials are evaluated with **Stirling's approximation** $n!\approx\sqrt{2\pi n}\,(n/e)^n$, which is already within $1\,\%$ for $n=10$.

> [!example] Worked Example 19.3 — Two dice
> Two fair dice are rolled. Find the probability that the sum is 7, and the probability that at least one die shows a 6.
>
> **Step 1 — sample space.** Ordered pairs $(d_1,d_2)$, $\lvert\Omega\rvert=6^2=36$, all equally likely. (Using *unordered* pairs, 21 outcomes, would be wrong: $\{3,4\}$ can occur in two ways but $\{3,3\}$ in only one, so unordered outcomes are not equally likely.)
>
> **Step 2 — sum 7.** $A=\{(1,6),(2,5),(3,4),(4,3),(5,2),(6,1)\}$, so $P(A)=6/36=1/6$.
>
> **Step 3 — at least one 6.** The complement is "no 6", which has $5\times5=25$ outcomes. Hence $P=1-25/36=11/36\approx0.306$. Inclusion–exclusion agrees: $1/6+1/6-1/36=11/36$. ✓

> [!example] Worked Example 19.4 — Sampling inspection (hypergeometric counting)
> A reel of 50 integrated circuits contains 5 defective ones. An inspector draws 4 at random without replacement. Find the probability that (a) none is defective, (b) exactly one is defective, (c) at least one is defective.
>
> **Step 1 — sample space.** Unordered samples of 4 from 50: $\binom{50}{4}=230\,300$, all equally likely.
>
> **Step 2 — none defective.** All 4 must come from the 45 good parts: $\binom{45}{4}=148\,995$. So
> $$
> P(\text{none})=\frac{\binom{45}{4}}{\binom{50}{4}}=\frac{148\,995}{230\,300}\approx0.6470 .
> $$
>
> **Step 3 — exactly one.** Choose 1 of the 5 defectives and 3 of the 45 good parts (product principle):
> $$
> P(\text{one})=\frac{\binom{5}{1}\binom{45}{3}}{\binom{50}{4}}=\frac{5\times14\,190}{230\,300}\approx0.3081 .
> $$
>
> **Step 4 — at least one.** $1-0.6470=0.3530$.
>
> In general, drawing $n$ from a population of $N$ containing $K$ "marked" items, the probability of exactly $k$ marked items is $\binom{K}{k}\binom{N-K}{n-k}/\binom{N}{n}$ — the **hypergeometric** distribution. Sampling inspection plans in quality control are designed with it ([Ch. 20](20-statistical-analysis.md)).

> [!example] Worked Example 19.5 — Address collisions (the birthday problem)
> Each of $n$ wireless sensor nodes picks a 16-bit address uniformly at random ($M=2^{16}=65\,536$ possibilities), independently of the others. What is the probability that at least two nodes pick the same address, for $n=100$ and $n=300$?
>
> **Step 1 — complement.** "No collision" means the $n$ addresses are all different. The number of ordered assignments with no repeat is $M(M-1)\cdots(M-n+1)$ out of $M^n$ equally likely assignments. So
> $$
> P(\text{no collision})=\prod_{k=0}^{n-1}\left(1-\frac{k}{M}\right).
> $$
>
> **Step 2 — approximation.** For $k\ll M$, $1-k/M\approx e^{-k/M}$, so the product is approximately $\exp\!\left(-\sum_{k=0}^{n-1}k/M\right)=\exp\!\left(-\dfrac{n(n-1)}{2M}\right)$.
>
> **Step 3 — numbers.** For $n=100$: exact $P(\text{collision})=0.0728$, approximation $1-e^{-9900/131072}=0.0727$. For $n=300$: exact $0.4961$, approximation $0.4956$.
>
> With only $300$ nodes — less than half a percent of the address space — a collision is a coin flip. The same calculation with $M=365$ gives the famous result that 23 people share a birthday with probability $0.507$. Collisions become likely once $n\approx\sqrt{M}$, which is why hash functions and random identifiers need about twice as many bits as one would naively expect.

> [!warning] Common Pitfall — "equally likely" must be justified
> The formula $P(A)=\lvert A\rvert/\lvert\Omega\rvert$ is valid only when the outcomes really are equally likely. Choosing a description in which they are not (unordered dice, "number of errors" as the outcome of a 3-bit word, "sum of the dice" as the outcome) and then counting gives wrong answers. When in doubt, describe outcomes at the finest level the physics makes symmetric — ordered tuples of individual trials — and count there.

> [!tip] Engineering Note — modelling first
> Almost every probability error in engineering reports is a modelling error, not an arithmetic one. Before computing, write down (1) the experiment, (2) the sample space, (3) the events of interest as subsets, (4) the source of the probability numbers (symmetry, data sheet, measurement). The axioms then do the rest mechanically.


## 19.2 Conditional probability

### Motivation

Information changes probabilities. The probability that a bit was received in error is $10^{-3}$ on average, but if the automatic gain control has just flagged a deep fade it may be $10^{-1}$. The probability that a power module fails in the next month is small, but much larger if it has already run hot for a year. Engineering decisions are almost always made *after* some partial information has arrived — a sensor reading, a test result, a received symbol — and conditional probability is the tool that updates the model accordingly.

### Definition

Suppose we learn that event $B$ has occurred. Outcomes outside $B$ are now impossible, and the outcomes inside $B$ should keep their *relative* likelihoods. The natural rule is therefore to restrict to $B$ and rescale so that $B$ has total probability 1.

> [!definition] Conditional probability
> For events $A,B$ with $P(B)>0$, the **conditional probability of $A$ given $B$** is
> $$
> P(A\mid B)=\frac{P(A\cap B)}{P(B)} .
> $$

In the relative-frequency picture: among the $N_B$ trials in which $B$ occurred, the fraction in which $A$ also occurred is $N_{AB}/N_B=(N_{AB}/N)/(N_B/N)\to P(AB)/P(B)$. For equally likely outcomes, $P(A\mid B)=\lvert A\cap B\rvert/\lvert B\rvert$: count inside the **reduced sample space** $B$.

> [!theorem] Conditional probability is a probability
> For fixed $B$ with $P(B)>0$, the map $A\mapsto P(A\mid B)$ satisfies all three Kolmogorov axioms on $(\Omega,\mathcal{F})$. Consequently every rule of Section 19.1 holds conditionally: $P(A^c\mid B)=1-P(A\mid B)$, $P(A_1\cup A_2\mid B)=P(A_1\mid B)+P(A_2\mid B)-P(A_1A_2\mid B)$, and so on.

*Proof.* Non-negativity is clear. $P(\Omega\mid B)=P(B)/P(B)=1$. If the $A_i$ are disjoint, so are the $A_i\cap B$, and $P\left(\bigcup A_i\mid B\right)=P\left(\bigcup(A_i\cap B)\right)/P(B)=\sum P(A_i\cap B)/P(B)=\sum P(A_i\mid B)$. $\square$

> [!example] Worked Example 19.6 — Reduced sample space
> Two fair dice are rolled. Given that the first die shows 6, what is the probability that the sum is at least 10? Compare with the unconditional probability.
>
> **Step 1 — condition.** $B=\{(6,1),\dots,(6,6)\}$, six equally likely outcomes.
>
> **Step 2 — count inside $B$.** The sum is at least 10 for $(6,4),(6,5),(6,6)$: $P(A\mid B)=3/6=1/2$.
>
> **Step 3 — unconditional.** Sums 10, 11, 12 occur in $3+2+1=6$ of 36 outcomes: $P(A)=1/6$.
>
> Knowing the first die is a 6 triples the probability. Check with the formula: $P(A\cap B)=3/36$, $P(B)=6/36$, ratio $1/2$. ✓

### The multiplication rule and the chain rule

Rearranging the definition gives the **multiplication rule**

$$
P(A\cap B)=P(A\mid B)\,P(B)=P(B\mid A)\,P(A),
$$

which is how joint probabilities are usually *computed* in sequential experiments: the probability of the first stage times the conditional probability of the second stage given the first. Applying it repeatedly yields the chain rule.

> [!theorem] Chain rule
> If $P(A_1A_2\cdots A_{n-1})>0$,
> $$
> P(A_1A_2\cdots A_n)=P(A_1)\,P(A_2\mid A_1)\,P(A_3\mid A_1A_2)\cdots P(A_n\mid A_1\cdots A_{n-1}).
> $$

*Proof.* Write each conditional probability as a ratio; the product telescopes: $P(A_1)\cdot\frac{P(A_1A_2)}{P(A_1)}\cdot\frac{P(A_1A_2A_3)}{P(A_1A_2)}\cdots=P(A_1\cdots A_n)$. $\square$

Sequential experiments are conveniently drawn as **probability trees**: each branch is labelled with a conditional probability given the path so far, and the probability of a complete path is the product of its branch labels. The labels on the branches leaving any node sum to 1, and the probabilities of all complete paths sum to 1 — two checks worth making every time.

> [!example] Worked Example 19.7 — Drawing without replacement (chain rule)
> A tray holds 10 MOSFETs, 3 of which are counterfeit. Parts are drawn one at a time without replacement. Find (a) the probability that the first two drawn are both counterfeit, and (b) the probability that the first counterfeit appears on the third draw.
>
> **(a)** Let $C_i$ = "draw $i$ is counterfeit". Then
> $$
> P(C_1C_2)=P(C_1)P(C_2\mid C_1)=\frac{3}{10}\cdot\frac{2}{9}=\frac{1}{15}\approx0.0667 .
> $$
> Counting check: $\binom{3}{2}/\binom{10}{2}=3/45=1/15$. ✓
>
> **(b)** The event is $C_1^cC_2^cC_3$. By the chain rule,
> $$
> P(C_1^cC_2^cC_3)=\frac{7}{10}\cdot\frac{6}{9}\cdot\frac{3}{8}=\frac{126}{720}=0.175 .
> $$
> Each factor is a conditional probability computed in the *current* tray: after two genuine parts are removed, 8 parts remain, 3 of them counterfeit.

### The law of total probability

Often $P(A)$ is hard to find directly but easy once we know which of several mutually exclusive "situations" holds. The law of total probability averages the conditional answers, weighting each situation by its probability.

> [!theorem] Law of total probability
> If $B_1,\dots,B_n$ is a partition of $\Omega$ with $P(B_i)>0$ for all $i$, then for any event $A$
> $$
> P(A)=\sum_{i=1}^{n}P(A\mid B_i)\,P(B_i).
> $$
> The same holds for a countably infinite partition.

*Proof.* Because the $B_i$ partition $\Omega$, the sets $A\cap B_i$ are disjoint and their union is $A$. By additivity $P(A)=\sum_iP(A\cap B_i)$, and each term equals $P(A\mid B_i)P(B_i)$ by the multiplication rule. $\square$

The simplest and most used case is the two-set partition $\{B,B^c\}$:

$$
P(A)=P(A\mid B)P(B)+P(A\mid B^c)\bigl(1-P(B)\bigr).
$$

Geometrically, $P(A)$ is a weighted average of the conditional probabilities $P(A\mid B_i)$, so it always lies between the smallest and the largest of them.

> [!example] Worked Example 19.8 — Three production lines
> A plant makes the same voltage regulator on three lines. Line 1 produces $50\,\%$ of the output with a defect rate of $1\,\%$, line 2 produces $30\,\%$ with $2\,\%$ defective, and line 3 produces $20\,\%$ with $3\,\%$ defective. What fraction of all regulators is defective?
>
> **Step 1 — partition.** $L_1,L_2,L_3$ (which line made the part), with $P(L_1)=0.5$, $P(L_2)=0.3$, $P(L_3)=0.2$; these sum to 1. ✓
>
> **Step 2 — conditionals.** $P(D\mid L_1)=0.01$, $P(D\mid L_2)=0.02$, $P(D\mid L_3)=0.03$.
>
> **Step 3 — total probability.**
> $$
> P(D)=0.01(0.5)+0.02(0.3)+0.03(0.2)=0.005+0.006+0.006=0.017 .
> $$
>
> So $1.7\,\%$ of regulators are defective — between the best line's $1\,\%$ and the worst's $3\,\%$, as it must be. Section 19.3 asks the reverse question: given a defective part, which line made it?

> [!example] Worked Example 19.9 — A channel with fades (Gilbert–Elliott model)
> A radio link is in a "good" state $90\,\%$ of the time, with bit-error probability $10^{-3}$, and in a "bad" (faded) state $10\,\%$ of the time, with bit-error probability $10^{-1}$. Find the average bit-error probability, and the probability that two consecutive bits (sent within the same state, the state lasting much longer than a bit) are both in error. Assume errors are independent *given the state*.
>
> **Step 1 — average error probability.** With $E$ = "bit in error", $G$ and $G^c$ the states:
> $$
> P(E)=10^{-3}(0.9)+10^{-1}(0.1)=0.0009+0.01=0.0109 .
> $$
> The bad state, present only $10\,\%$ of the time, contributes over $90\,\%$ of the errors.
>
> **Step 2 — two errors.** Conditionally on the state, two errors have probability $p^2$:
> $$
> P(E_1E_2)=(10^{-3})^2(0.9)+(10^{-1})^2(0.1)=9\times10^{-7}+10^{-3}\approx1.0009\times10^{-3}.
> $$
>
> **Step 3 — compare with a memoryless channel.** If bits were independent with the average error rate, $P(E_1E_2)$ would be $0.0109^2=1.19\times10^{-4}$ — more than eight times too small. Moreover $P(E_2\mid E_1)=1.0009\times10^{-3}/0.0109\approx0.092$: once one error is seen, the next bit is almost 10 times more likely than average to be wrong. Errors on fading channels come in **bursts**, which is why such links use interleaving before the error-correcting code.

> [!warning] Common Pitfall — the conditioning bar is not symmetric
> $P(A\mid B)$ and $P(B\mid A)$ are different numbers, related by $P(A\mid B)P(B)=P(B\mid A)P(A)$. "Most failed units had overheated" ($P(\text{hot}\mid\text{fail})$ large) does not mean "most units that overheat fail" ($P(\text{fail}\mid\text{hot})$ may be small). Confusing the two — sometimes called the *prosecutor's fallacy* — is the most common reasoning error in fault analysis. Bayes' theorem (Section 19.3) converts one into the other.

> [!warning] Common Pitfall — averaging conditionals with the wrong weights
> The law of total probability weights $P(A\mid B_i)$ by $P(B_i)$, not equally. Averaging the three defect rates in Example 19.8 without weights gives $2\,\%$, not the correct $1.7\,\%$. And the partition must be genuinely exhaustive and disjoint; overlapping "situations" double-count.

> [!tip] Engineering Note — conditioning as a design tool
> The law of total probability is the formal version of "divide into cases". It appears as the averaging of error rates over fading states, of failure rates over operating temperatures, of delay over network loads. When a direct calculation looks hard, ask: *which unknown, if I knew it, would make the problem easy?* Condition on it, then average.

## 19.3 Bayes' theorem

### Motivation

Physical models usually run from cause to effect. A data sheet tells us the probability that a detector raises an alarm *given* that a fault is present; a channel model tells us the probability of receiving a 1 *given* that a 0 was sent. But the engineer's question runs the other way: the alarm has sounded — is there a fault? A 1 has been received — what was sent? Bayes' theorem inverts the conditioning, turning cause-to-effect models into effect-to-cause inferences. It is the mathematical core of detection, diagnosis, classification and estimation.

### Statement and derivation

> [!theorem] Bayes' theorem
> Let $B_1,\dots,B_n$ be a partition of $\Omega$ with $P(B_i)>0$, and let $A$ be an event with $P(A)>0$. Then
> $$
> P(B_k\mid A)=\frac{P(A\mid B_k)\,P(B_k)}{\displaystyle\sum_{i=1}^{n}P(A\mid B_i)\,P(B_i)} .
> $$
> For a two-set partition $\{B,B^c\}$:
> $$
> P(B\mid A)=\frac{P(A\mid B)P(B)}{P(A\mid B)P(B)+P(A\mid B^c)P(B^c)} .
> $$

*Derivation.* The multiplication rule gives $P(B_k\cap A)$ in two ways: $P(B_k\mid A)P(A)=P(A\mid B_k)P(B_k)$. Divide by $P(A)$ and expand $P(A)$ with the law of total probability. $\square$

The terms have standard names, which carry over unchanged to estimation theory ([Ch. 20](20-statistical-analysis.md)):

| Term | Name | Meaning |
|---|---|---|
| $P(B_k)$ | **prior** | belief in hypothesis $B_k$ before the observation |
| $P(A\mid B_k)$ | **likelihood** | how probable the observation is if $B_k$ is true |
| $P(A)$ | **evidence** (normalizer) | total probability of the observation |
| $P(B_k\mid A)$ | **posterior** | updated belief after observing $A$ |

In words: *posterior is proportional to likelihood times prior*, and the denominator merely makes the posteriors sum to 1. This suggests a practical recipe: compute the products $P(A\mid B_i)P(B_i)$ for every hypothesis, then divide each by their sum.

> [!example] Worked Example 19.10 — Which line made the defective part?
> Continue Example 19.8. A regulator selected at random is found defective. Find the probability that it came from each line.
>
> **Step 1 — joint products.** $P(D\mid L_i)P(L_i)$: line 1, $0.005$; line 2, $0.006$; line 3, $0.006$. Sum $=P(D)=0.017$.
>
> **Step 2 — normalize.**
> $$
> P(L_1\mid D)=\frac{0.005}{0.017}=0.294,\qquad P(L_2\mid D)=\frac{0.006}{0.017}=0.353,\qquad P(L_3\mid D)=\frac{0.006}{0.017}=0.353 .
> $$
> Check: the posteriors sum to 1. ✓
>
> Line 1 makes half the parts but only $29\,\%$ of the defectives; line 3 makes a fifth of the parts but $35\,\%$ of the defectives. The prior (production share) and the likelihood (defect rate) pull in opposite directions, and Bayes' theorem weighs them exactly.

### The binary channel

The simplest communication model has a binary input $X\in\{0,1\}$ and a binary output $Y\in\{0,1\}$. It is specified by the **prior** probabilities $P(X=0)=p_0$, $P(X=1)=p_1=1-p_0$, and the **transition probabilities** $P(Y=y\mid X=x)$. In the **binary symmetric channel** (BSC) both kinds of error have the same probability $\varepsilon$: $P(Y=1\mid X=0)=P(Y=0\mid X=1)=\varepsilon$. In an asymmetric channel (optical on–off keying, where a "1" can be lost but noise rarely creates one, is a typical example) the two error probabilities differ.

> [!example] Worked Example 19.11 — Binary asymmetric channel
> A source sends $0$ with probability $0.6$ and $1$ with probability $0.4$. The channel has $P(Y=1\mid X=0)=0.05$ and $P(Y=0\mid X=1)=0.10$. Find (a) the output probabilities, (b) the posterior probabilities of the input given each output, (c) the error probability if the receiver simply decides $\hat X=Y$.
>
> **(a) Total probability.**
> $$
> P(Y=1)=0.05(0.6)+0.90(0.4)=0.03+0.36=0.39,\qquad P(Y=0)=0.95(0.6)+0.10(0.4)=0.57+0.04=0.61 .
> $$
>
> **(b) Bayes.**
> $$
> P(X=1\mid Y=1)=\frac{0.36}{0.39}=0.923,\qquad P(X=0\mid Y=0)=\frac{0.57}{0.61}=0.934 .
> $$
> So a received 1 is wrong with probability $0.077$ and a received 0 with probability $0.066$.
>
> **(c) Error probability.** An error occurs when $(X,Y)=(0,1)$ or $(1,0)$:
> $$
> P_e=P(Y=1\mid X=0)p_0+P(Y=0\mid X=1)p_1=0.05(0.6)+0.10(0.4)=0.07 .
> $$
> Check: $P_e=P(Y=1)(1-0.923)+P(Y=0)(1-0.934)=0.39(0.0769)+0.61(0.0656)=0.030+0.040=0.070$. ✓

### MAP and ML decisions

A receiver that observes $Y=y$ must decide which input was sent. The decision that minimizes the probability of error is the one that picks the input with the **largest posterior**.

> [!theorem] MAP decision rule
> Given an observation $A$ and hypotheses $B_1,\dots,B_n$, the rule "decide $B_k$ where $k$ maximizes $P(B_k\mid A)$" — equivalently, maximizes $P(A\mid B_k)P(B_k)$ — minimizes the probability of a wrong decision. It is called the **maximum a posteriori (MAP)** rule. When all priors are equal, it reduces to the **maximum-likelihood (ML)** rule: maximize $P(A\mid B_k)$.

*Proof sketch.* For any decision rule, the probability of a *correct* decision is $\sum_{\text{observations } a}P(\text{decided hypothesis is true}\mid a)P(a)$. Each term is maximized separately by choosing, for each $a$, the hypothesis with the largest posterior. Since the denominators $P(a)$ are common to all hypotheses, comparing posteriors is the same as comparing the products likelihood $\times$ prior. $\square$

For two hypotheses the MAP rule is a **likelihood-ratio test**:

$$
\Lambda(a)=\frac{P(a\mid B_1)}{P(a\mid B_0)}\ \underset{B_0}{\overset{B_1}{\gtrless}}\ \frac{P(B_0)}{P(B_1)} ,
$$

i.e. decide $B_1$ if the likelihood ratio exceeds the prior odds against it. Section 19.7 applies exactly this test to a signal in Gaussian noise.

> [!example] Worked Example 19.12 — When MAP ignores the channel
> In the channel of Example 19.11, MAP gives: for $Y=0$, compare $0.95(0.6)=0.57$ with $0.10(0.4)=0.04$ — decide 0; for $Y=1$, compare $0.05(0.6)=0.03$ with $0.90(0.4)=0.36$ — decide 1. So $\hat X=Y$ is already MAP, and $P_e=0.07$ is the minimum possible.
>
> Now suppose 1s are rare: $p_1=0.03$, $p_0=0.97$, same channel.
>
> **Step 1 — $Y=1$ received.** Compare $P(Y=1\mid 0)p_0=0.05(0.97)=0.0485$ with $P(Y=1\mid 1)p_1=0.90(0.03)=0.027$. The first is larger: **decide 0 even though 1 was received.**
>
> **Step 2 — $Y=0$ received.** $0.95(0.97)$ vastly exceeds $0.10(0.03)$: decide 0.
>
> **Step 3 — error probabilities.** The MAP receiver always outputs 0 and is wrong only when a 1 was sent: $P_e=0.03$. The "trust the channel" receiver $\hat X=Y$ has $P_e=0.05(0.97)+0.10(0.03)=0.0515$.
>
> With strongly unequal priors, a noisy observation may carry less information than the prior. In practice one would not build such a receiver — it never detects a 1 — which shows that "minimum error probability" is not always the right criterion. Detection theory therefore also uses rules that weight the two errors by different *costs* (Bayes risk) or fix the false-alarm rate (Neyman–Pearson).

### Base rates and fault detection

Bayes' theorem explains why detectors for rare events produce mostly false alarms, a fact that surprises almost everyone the first time.

> [!example] Worked Example 19.13 — Insulation fault detector
> An online partial-discharge monitor is fitted to a population of cable joints, of which $0.2\,\%$ have an incipient insulation fault. The monitor detects a real fault with probability $0.99$ (sensitivity) and raises a false alarm on a healthy joint with probability $0.02$. (a) An alarm is raised on a joint. What is the probability that it is faulty? (b) A second, independent monitor of the same specification also alarms on that joint. What now?
>
> **Step 1 — set up.** $F$ = faulty, $A$ = alarm. Prior $P(F)=0.002$. Likelihoods $P(A\mid F)=0.99$, $P(A\mid F^c)=0.02$.
>
> **Step 2 — Bayes.**
> $$
> P(F\mid A)=\frac{0.99(0.002)}{0.99(0.002)+0.02(0.998)}=\frac{0.00198}{0.00198+0.01996}=0.0902 .
> $$
> Despite a "99 % accurate" detector, only about **9 %** of alarms are real. Out of 100 000 joints, 200 are faulty and 198 of those alarm; but 2% of the 99 800 healthy joints — 1996 — also alarm.
>
> **Step 3 — odds form for the second alarm.** Write Bayes in odds: posterior odds = likelihood ratio $\times$ prior odds, where the likelihood ratio of an alarm is $\Lambda=0.99/0.02=49.5$. After the first alarm the odds are $0.0902/0.9098=0.0992$. A second alarm, conditionally independent given the joint's state, multiplies the odds again by 49.5:
> $$
> \text{odds}=0.0992\times49.5=4.91,\qquad P(F\mid A_1A_2)=\frac{4.91}{1+4.91}=0.831 .
> $$
>
> Two independent alarms move the probability from $0.2\,\%$ to $83\,\%$. Confirmation by an independent test is the standard engineering cure for the base-rate problem.

The **odds form** used in Step 3 is worth stating separately because it makes sequential updating trivial:

$$
\underbrace{\frac{P(B\mid A)}{P(B^c\mid A)}}_{\text{posterior odds}}=\underbrace{\frac{P(A\mid B)}{P(A\mid B^c)}}_{\text{likelihood ratio}}\times\underbrace{\frac{P(B)}{P(B^c)}}_{\text{prior odds}} .
$$

With independent observations the likelihood ratios multiply; taking logarithms, the evidence *adds*, which is why practical detectors accumulate **log-likelihood ratios**.

The dependence of the posterior on the prior is shown below for the detector of Example 19.13.

```desmos-graph
left=0; right=0.2; bottom=0; top=1.05
---
y=\frac{0.99x}{0.99x+0.02(1-x)}|blue
y=\frac{0.99^{2}x}{0.99^{2}x+0.02^{2}(1-x)}|red|dashed
(0.002,0.0902)|label:prior 0.002
(0.01,0.3333)|label:prior 0.01
```

*Posterior probability of a fault after one alarm (blue) and after two independent alarms (red, dashed) as a function of the prior fault probability $x$, for sensitivity 0.99 and false-alarm rate 0.02. For rare faults the single-alarm posterior collapses towards zero; a second alarm restores confidence.*

> [!warning] Common Pitfall — base-rate neglect
> A test's sensitivity $P(A\mid F)$ and specificity $1-P(A\mid F^c)$ describe the *test*; the probability that an alarm is real, $P(F\mid A)$, depends also on the *population* through the prior $P(F)$. Quoting "99 % accurate" without the base rate says almost nothing about how alarms should be treated. The same detector deployed on a population with $10\,\%$ faulty units has $P(F\mid A)=0.85$.

> [!tip] Engineering Note — Bayes everywhere
> The binary channel and the fault detector are the same mathematics. So are spam filters, radar target detection, built-in self-test diagnostics, and the decoding of error-correcting codes (where the "hypotheses" are the possible codewords). In every case: write the prior, write the likelihood of the observation under each hypothesis, multiply, normalize, and decide by the largest posterior or by a cost-weighted threshold.

## 19.4 Independence

### Motivation

Most probability models of real systems are built from parts: the bits of a word, the components of a circuit, the successive packets on a link, the samples of a noise waveform. The model becomes tractable when the parts do not influence each other, because then the probability of a compound event is simply the product of the probabilities of its pieces. Independence is that "no influence" assumption made precise. It is the most powerful simplification in applied probability — and the most frequently abused one, so knowing exactly what it means, and when it fails, matters as much as using it.

### Independence of two events

Intuitively, $A$ and $B$ are independent if learning that $B$ occurred does not change the probability of $A$: $P(A\mid B)=P(A)$. Multiplying by $P(B)$ gives a symmetric form that also makes sense when $P(B)=0$, and this is taken as the definition.

> [!definition] Independent events
> Events $A$ and $B$ are **independent** if
> $$
> P(A\cap B)=P(A)\,P(B).
> $$
> If $P(B)>0$ this is equivalent to $P(A\mid B)=P(A)$; if $P(A)>0$, to $P(B\mid A)=P(B)$.

> [!theorem] Complements of independent events
> If $A$ and $B$ are independent, so are $A$ and $B^c$, $A^c$ and $B$, and $A^c$ and $B^c$.

*Proof.* $P(A\cap B^c)=P(A)-P(A\cap B)=P(A)-P(A)P(B)=P(A)\bigl(1-P(B)\bigr)=P(A)P(B^c)$. The other cases follow by applying this twice. $\square$

> [!warning] Common Pitfall — independent is not the same as disjoint
> Disjoint events with non-zero probabilities are never independent: if $A\cap B=\varnothing$ then $P(AB)=0\ne P(A)P(B)$. Indeed, disjoint events are strongly *dependent* — knowing that $A$ occurred tells you for certain that $B$ did not. "Mutually exclusive" is a property of the sets; "independent" is a property of the probability measure.

### Independence of several events

For three or more events, pairwise independence is not enough.

> [!definition] Mutual independence
> Events $A_1,\dots,A_n$ are **(mutually) independent** if for *every* subcollection $A_{i_1},\dots,A_{i_k}$ with $2\le k\le n$ and distinct indices,
> $$
> P(A_{i_1}\cap A_{i_2}\cap\cdots\cap A_{i_k})=P(A_{i_1})P(A_{i_2})\cdots P(A_{i_k}).
> $$
> For $n$ events this is $2^n-n-1$ conditions. An infinite family is independent if every finite subfamily is.

> [!example] Worked Example 19.14 — Pairwise but not mutually independent (parity bit)
> Two fair, independent data bits $b_1,b_2$ are sent with an even-parity bit $b_3=b_1\oplus b_2$. Let $A_i$ = "$b_i=1$". Are $A_1,A_2,A_3$ independent?
>
> **Step 1 — marginals.** The four equally likely words are $000,011,101,110$. Each bit equals 1 in exactly two of them: $P(A_i)=1/2$.
>
> **Step 2 — pairs.** $A_1\cap A_2=\{110\}$, $A_1\cap A_3=\{101\}$, $A_2\cap A_3=\{011\}$, each with probability $1/4=P(A_i)P(A_j)$. So the events are **pairwise independent**.
>
> **Step 3 — triple.** $A_1\cap A_2\cap A_3=\varnothing$ (no valid word is 111), so $P=0\ne1/8$. The events are **not mutually independent**: any two bits determine the third.
>
> This is exactly why a parity bit detects errors: the three bits are not free.

### Conditional independence

Events $A$ and $B$ are **conditionally independent given $C$** if $P(AB\mid C)=P(A\mid C)P(B\mid C)$. Conditional independence neither implies nor is implied by independence. The fading channel of Example 19.9 is the standard illustration: given the channel state, bit errors are independent; *unconditionally* they are strongly dependent, because an error is evidence that the channel is in the bad state, which makes the next error more likely. The second alarm in Example 19.13 was assumed conditionally independent given the joint's state — not independent outright, since both alarms respond to the same fault.

### Repeated independent trials

The most important use of independence is to build the probability law of $n$ repeated, independent trials of the same experiment. Suppose each trial results in "success" with probability $p$ and "failure" with probability $q=1-p$ (a **Bernoulli trial**). By independence, any particular sequence with $k$ successes and $n-k$ failures has probability $p^kq^{n-k}$, and there are $\binom{n}{k}$ such sequences (choose the positions of the successes). Therefore:

> [!theorem] Binomial probability law
> In $n$ independent Bernoulli trials with success probability $p$, the probability of exactly $k$ successes is
> $$
> P(k\text{ successes})=\binom{n}{k}p^k(1-p)^{n-k},\qquad k=0,1,\dots,n .
> $$
> These probabilities sum to $(p+q)^n=1$ by the binomial theorem.

> [!example] Worked Example 19.15 — Bit errors and a Hamming code
> A binary channel has independent bit errors with probability $p=0.01$. (a) What is the probability that an uncoded 4-bit data word is received correctly? (b) The 4 data bits are instead sent as a 7-bit Hamming codeword, which corrects any single error. What is the probability that the decoded word is wrong?
>
> **(a)** All 4 bits correct: $(0.99)^4=0.9606$; word error probability $1-0.9606=0.0394$.
>
> **(b)** The decoder fails only if 2 or more of the 7 bits are wrong:
> $$
> P(\text{fail})=1-\left[\binom{7}{0}(0.99)^7+\binom{7}{1}(0.01)(0.99)^6\right]=1-[0.93207+0.06590]=0.00203 .
> $$
>
> The code reduces the word-error probability by a factor of about 19, at the price of sending 7 bits instead of 4. The dominant failure term is $\binom{7}{2}p^2q^5\approx 0.00199$, so $P(\text{fail})\approx 21p^2$: a single-error-correcting code turns an error rate proportional to $p$ into one proportional to $p^2$. This scaling argument is the essence of coding gain.

### Reliability of systems

**Reliability** $R$ is the probability that a component or system performs its function (over a stated time under stated conditions). If components fail independently, the reliability of a system follows from its structure by the rules of this section. Let component $i$ work with probability $p_i$ and fail with probability $q_i=1-p_i$.

**Series system.** The system works only if *every* component works (a chain: fuse, converter, regulator, load). By independence,

$$
R_{\text{series}}=p_1p_2\cdots p_n=\prod_{i=1}^{n}p_i .
$$

A series system is less reliable than its weakest component, and reliability decays geometrically with the number of parts: 100 parts of reliability $0.999$ give only $0.905$.

**Parallel (redundant) system.** The system works if *at least one* component works. By De Morgan, it fails only if all fail:

$$
R_{\text{parallel}}=1-q_1q_2\cdots q_n=1-\prod_{i=1}^{n}(1-p_i) .
$$

A parallel system is more reliable than its best component. Two units of $0.9$ in parallel give $0.99$; three give $0.999$ — each additional unit adds "one nine".

**$k$-out-of-$n$ system.** The system works if at least $k$ of $n$ identical, independent components work. By the binomial law,

$$
R_{k/n}=\sum_{i=k}^{n}\binom{n}{i}p^i(1-p)^{n-i} .
$$

Series is $n$-out-of-$n$, parallel is 1-out-of-$n$. The important special case 2-out-of-3 is **triple modular redundancy (TMR)** with a majority voter: $R_{\text{TMR}}=3p^2-2p^3$ (for a perfect voter).

**Series–parallel networks** are reduced step by step: replace each parallel group by an equivalent component with the parallel formula, each series chain by the series formula, and repeat — much as series and parallel resistors are combined. Networks that cannot be reduced this way (the **bridge** is the classic one) are handled by **conditioning on a key component** with the law of total probability:

$$
R=R(\text{system}\mid\text{component } k\text{ works})\,p_k+R(\text{system}\mid\text{component } k\text{ fails})\,q_k ,
$$

where in each case the network simplifies: a working component is replaced by a short circuit, a failed one by an open circuit.

> [!example] Worked Example 19.16 — A redundant power path
> A board is powered through an input fuse ($p=0.999$), then two identical DC–DC converters in parallel (each $p=0.95$), then a supervisory controller ($p=0.99$). Find the system reliability, and compare with a single converter.
>
> **Step 1 — reduce the parallel pair.** $R_{\text{conv}}=1-(0.05)^2=0.9975$.
>
> **Step 2 — series chain.** $R=0.999\times0.9975\times0.99=0.98654$.
>
> **Step 3 — without redundancy.** $R=0.999\times0.95\times0.99=0.93955$.
>
> The second converter cuts the failure probability from $6.0\,\%$ to $1.35\,\%$. The fuse and controller now dominate: the remaining unreliability is about $0.1\,\%+0.25\,\%+1\,\%$. Further redundancy should go where the largest series term is — here the controller.

> [!example] Worked Example 19.17 — 2-out-of-3 voting
> Three channels of a safety sensor each work with probability $0.95$, independently; a perfect voter outputs the majority. Find the system reliability.
>
> **Step 1 — binomial.** At least 2 of 3 must work:
> $$
> R=\binom{3}{2}(0.95)^2(0.05)+\binom{3}{3}(0.95)^3=0.135375+0.857375=0.99275 .
> $$
>
> **Step 2 — check with the formula.** $3p^2-2p^3=3(0.9025)-2(0.857375)=2.7075-1.71475=0.99275$. ✓
>
> TMR improves on a single channel ($0.95$) — but only because $p>1/2$. Solving $3p^2-2p^3=p$ gives $p=1/2$: if each channel is more likely to fail than to work, majority voting makes things *worse*.

> [!example] Worked Example 19.18 — The bridge network
> Five links connect a source $S$ to a terminal $T$: link 1 from $S$ to node $u$, link 2 from $u$ to $T$, link 3 from $S$ to node $v$, link 4 from $v$ to $T$, and a bridging link 5 between $u$ and $v$. Each works independently with probability $p=0.9$. Find the probability that $S$ and $T$ are connected.
>
> **Step 1 — condition on link 5 working.** Nodes $u$ and $v$ merge. The network becomes ($1\parallel3$) in series with ($2\parallel4$):
> $$
> R_1=\bigl[1-(0.1)^2\bigr]^2=0.99^2=0.9801 .
> $$
>
> **Step 2 — condition on link 5 failed.** The network is two series paths in parallel, $(1\text{–}2)\parallel(3\text{–}4)$:
> $$
> R_0=1-(1-0.81)^2=1-0.0361=0.9639 .
> $$
>
> **Step 3 — total probability.**
> $$
> R=0.9(0.9801)+0.1(0.9639)=0.88209+0.09639=0.97848 .
> $$
>
> A brute-force check — summing the probabilities of all $2^5=32$ working/failed patterns for which a path exists — gives the same $0.97848$. Conditioning reduced 32 cases to 2.

### Reliability as a function of time

If a component's lifetime $T$ has the exponential distribution (Section 19.6) with constant failure rate $\lambda$, its reliability at time $t$ is $R(t)=P(T>t)=e^{-\lambda t}$. Substituting $p=e^{-\lambda t}$ into the structure formulas gives the reliability of a system as a function of time:

$$
R_{\text{series}}(t)=e^{-n\lambda t},\qquad R_{\text{par},2}(t)=2e^{-\lambda t}-e^{-2\lambda t},\qquad R_{\text{TMR}}(t)=3e^{-2\lambda t}-2e^{-3\lambda t}.
$$

A series system of exponential components is again exponential, with failure rate equal to the sum of the component rates — the basis of the "parts count" reliability prediction method.

```desmos-graph
left=0; right=3; bottom=0; top=1.05
---
y=e^{-x}|black
y=e^{-2x}|red|dashed
y=2e^{-x}-e^{-2x}|blue
y=3e^{-2x}-2e^{-3x}|green
(0.6931,0.5)|label:crossover at ln 2
```

*Reliability versus normalized time $\lambda t$ for a single unit (black), two units in series (red, dashed), two in parallel (blue) and 2-out-of-3 TMR (green). TMR is better than a single unit early in life but worse after $\lambda t=\ln 2$, when each channel's reliability drops below $1/2$; it suits short missions with repair or reset.*

> [!warning] Common Pitfall — common-cause failures
> The product formulas assume *independent* failures. Two "redundant" converters fed from the same overheated enclosure, two sensors calibrated with the same faulty reference, or three TMR channels running the same buggy software are not independent. A common cause that fails all channels at once with probability $10^{-3}$ caps the system reliability at about $0.999$ no matter how much redundancy is added. Reliability engineers model this explicitly (the $\beta$-factor model) or diversify the redundant channels.

> [!tip] Engineering Note — where independence comes from
> Independence is a *physical* assumption, justified by separation of mechanisms: different thermal noise sources in different resistors, bits separated by more than the channel's coherence time, components that share no stress. Write down why you believe it before using the product rule, and check it against data when you can — Ch. 20 gives tests for exactly this.

## 19.5 Random variables

### Motivation

Outcomes of experiments are often not numbers: a received word, the state of a network, a sequence of heads and tails. Engineering questions, however, are almost always numerical — *how many* bits are wrong, *what* voltage was measured, *how long* until failure, *how much* power is in the noise. A random variable is the device that attaches a number to each outcome. Once we have numbers we can add, average, threshold and filter them, and the full machinery of calculus becomes available. From here on, almost all of probability is phrased in terms of random variables, and the underlying sample space fades into the background.

### Definition

> [!definition] Random variable
> A **random variable** (RV) on a probability space $(\Omega,\mathcal{F},P)$ is a function $X:\Omega\to\mathbb{R}$ that assigns a real number $X(\zeta)$ to each outcome $\zeta$, such that for every real $x$ the set
> $$
> \{X\le x\}=\{\zeta\in\Omega: X(\zeta)\le x\}
> $$
> is an event (belongs to $\mathcal{F}$). The set of values $X$ can take is its **range** $S_X$.

The technical condition (called *measurability*) guarantees that we can ask for $P(X\le x)$; it is automatically satisfied by every function an engineer writes down, and we shall not mention it again. Note the language: a random variable is neither random nor a variable — it is a deterministic *function*. The randomness lives in which outcome $\zeta$ occurs.

Once $\{X\le x\}$ is an event for every $x$, so are all the sets built from them by countable unions, intersections and complements, for instance

$$
\{a<X\le b\}=\{X\le b\}\setminus\{X\le a\},\qquad \{X=a\}=\bigcap_{n=1}^{\infty}\{a-1/n<X\le a\},\qquad \{X>a\}=\{X\le a\}^c .
$$

So every question of the form "does $X$ fall in the set $B$?" has a probability $P(X\in B)$. The assignment $B\mapsto P(X\in B)$ is itself a probability measure on the real line, called the **distribution** (or **law**) of $X$. Two random variables defined on completely different sample spaces can have the same distribution; for most purposes they are then interchangeable.

> [!example] Worked Example 19.19 — Number of bit errors
> A 3-bit word is sent over a channel with independent bit errors of probability $p=0.1$. Let $X$ = number of bits in error. Find the probability of each value of $X$.
>
> **Step 1 — the function.** $\Omega=\{000,001,\dots,111\}$ (1 = error) and $X(\zeta)$ is the number of 1s in $\zeta$: $X(000)=0$, $X(010)=1$, $X(111)=3$, and so on. $S_X=\{0,1,2,3\}$.
>
> **Step 2 — outcome probabilities.** By independence, an outcome with $k$ ones has probability $0.1^k\,0.9^{3-k}$.
>
> **Step 3 — collect outcomes with the same value of $X$.**
> $$
> \begin{aligned}
> P(X=0)&=0.9^3=0.729, & P(X=1)&=3(0.1)(0.9)^2=0.243,\\
> P(X=2)&=3(0.1)^2(0.9)=0.027, & P(X=3)&=0.1^3=0.001 .
> \end{aligned}
> $$
> The probabilities sum to 1. ✓ Eight outcomes have been condensed into four numbers — the information relevant to the question "how many errors?".

### Types of random variable

> [!definition] Discrete, continuous and mixed random variables
> - $X$ is **discrete** if its range $S_X$ is finite or countably infinite: $X$ takes isolated values $x_1,x_2,\dots$ with probabilities $P(X=x_k)>0$ that sum to 1. (Counts, symbols, bits, quantizer outputs.)
> - $X$ is **continuous** if $P(X=x)=0$ for every single value $x$ and probabilities of intervals are given by integrals of a density: $P(a<X\le b)=\int_a^b f_X(x)\,dx$. (Voltages, times, phases, positions.)
> - $X$ is **mixed** if it has both: some isolated values with positive probability, and a continuous part spread over intervals. (A clipped voltage, the waiting time of a packet that is either sent at once or queued.)

Every random variable is completely described by its **cumulative distribution function** $F_X(x)=P(X\le x)$ (Section 19.8). For discrete RVs it is more convenient to list the probabilities of the individual values, the **probability mass function** (PMF); for continuous RVs, the derivative of the CDF, the **probability density function** (PDF). Both are developed in Section 19.7.

> [!warning] Common Pitfall — "probability zero" is not "impossible"
> For a continuous RV every individual value has probability zero: the thermal-noise voltage is exactly $1.000\ldots\,\text{mV}$ with probability 0, yet some value always occurs. Zero-probability events can happen; they are just too "thin" to carry probability. Consequently, for a continuous RV, $P(a<X<b)=P(a\le X\le b)$ — endpoints do not matter. For discrete and mixed RVs they matter a great deal.

### Indicator random variables

The simplest random variable takes only the values 0 and 1.

> [!definition] Indicator random variable
> For an event $A$, the **indicator** of $A$ is
> $$
> I_A(\zeta)=\begin{cases}1,&\zeta\in A\\0,&\zeta\notin A.\end{cases}
> $$
> It satisfies $P(I_A=1)=P(A)$, $I_{A^c}=1-I_A$, $I_{A\cap B}=I_AI_B$, and $I_{A\cup B}=I_A+I_B-I_AI_B$.

Indicators turn events into numbers and set operations into arithmetic. Their great use is in **counting**: the number of events among $A_1,\dots,A_n$ that occur is $X=I_{A_1}+\cdots+I_{A_n}$. In Example 19.19, $X$ is the sum of three indicators "bit $i$ is wrong". Section 19.9 shows that this decomposition makes the mean of $X$ trivial to compute even when the events are dependent.

### Functions of a random variable

If $X$ is a random variable and $g$ is an ordinary real function, then $Y=g(X)$ — the function $\zeta\mapsto g(X(\zeta))$ — is again a random variable. This is the probabilistic version of passing a signal through a memoryless system: a diode, a square-law detector, a limiter, a quantizer, an amplifier. The central question is: **given the distribution of $X$, what is the distribution of $Y$?**

The principle is always the same: an event about $Y$ is an event about $X$,

$$
P(Y\in B)=P\bigl(g(X)\in B\bigr)=P\bigl(X\in g^{-1}(B)\bigr),\qquad g^{-1}(B)=\{x: g(x)\in B\}.
$$

For a discrete $X$ this simply means collecting probability: every value $y$ of $Y$ gets the total probability of the $x$-values that map to it,

$$
P(Y=y)=\sum_{x:\ g(x)=y}P(X=x).
$$

For continuous $X$ the same principle leads to the CDF method and the density-transformation formula of Sections 19.7 and 19.8.

> [!example] Worked Example 19.20 — A squarer acting on a discrete signal
> A 5-level signal $X$ takes the values $-2,-1,0,1,2$ V with probabilities $0.1,0.2,0.4,0.2,0.1$. It drives a square-law power detector whose output is $Y=X^2$ (in V$^2$). Find the PMF of $Y$.
>
> **Step 1 — range of $Y$.** $g(x)=x^2$ maps $\{-2,-1,0,1,2\}$ onto $\{0,1,4\}$.
>
> **Step 2 — collect preimages.**
> $$
> P(Y=0)=P(X=0)=0.4,\qquad P(Y=1)=P(X=-1)+P(X=1)=0.4,\qquad P(Y=4)=P(X=-2)+P(X=2)=0.2 .
> $$
>
> **Step 3 — check.** $0.4+0.4+0.2=1$. ✓ Because $g$ is not one-to-one, information is lost: from $Y=1$ we cannot tell whether $X$ was $+1$ or $-1$.

> [!example] Worked Example 19.21 — A limiter produces a mixed random variable
> A voltage $X$ is equally likely to be anywhere in $(-2,2)$ V, i.e. $P(a<X<b)=(b-a)/4$ for $-2\le a<b\le 2$. It passes through an ideal limiter with output $Y=X$ for $\lvert X\rvert\le1$, $Y=+1$ for $X>1$ and $Y=-1$ for $X<-1$. Describe the distribution of $Y$.
>
> **Step 1 — clipped values.** $P(Y=1)=P(X>1)=(2-1)/4=0.25$; likewise $P(Y=-1)=0.25$. These are point masses.
>
> **Step 2 — the linear region.** For $-1<a<b<1$, $P(a<Y<b)=P(a<X<b)=(b-a)/4$: probability is spread uniformly over $(-1,1)$ with density $1/4$ per volt, total $0.5$.
>
> **Step 3 — check.** $0.25+0.25+0.5=1$. ✓
>
> $Y$ is **mixed**: it is neither discrete (it takes a continuum of values) nor continuous ($P(Y=1)>0$). Saturating amplifiers, ADCs driven beyond full scale, and rectifiers all produce mixed outputs; Section 19.8 draws their staircase-plus-ramp CDFs, and Section 19.7 represents the point masses with Dirac impulses.

> [!tip] Engineering Note — random variables are measurements
> A good habit is to read "random variable" as "the reading on an instrument attached to the experiment". Different instruments on the same experiment give different RVs on the same sample space ($X$ = number of errors, $Y$ = position of the first error, $Z$ = parity of the received word); their joint behaviour is the subject of the joint distributions introduced in Section 19.6. And a memoryless circuit applied to a reading gives a new reading, $Y=g(X)$, whose distribution is fixed by that of $X$.

## 19.6 Probability distributions

### Motivation

A handful of distributions account for the overwhelming majority of probability models in electrical engineering. Each arises from a simple mechanism: counting successes in independent trials (binomial), waiting for the first success (geometric), counting rare events in a long interval (Poisson), complete ignorance over an interval (uniform), memoryless lifetimes (exponential), the sum of many small effects (Gaussian), the envelope of a two-dimensional Gaussian (Rayleigh). Recognizing the mechanism lets you write down the distribution without collecting data. This section introduces these key families and the first ideas of **joint** distributions of two random variables; [Ch. 20](20-statistical-analysis.md) extends the catalogue (chi-square, Student $t$, lognormal, Weibull, …) and shows how to fit and test distributions against data.

> [!definition] Distribution of a random variable
> The **distribution** of $X$ is the rule that gives $P(X\in B)$ for every event $B\subseteq\mathbb{R}$. In practice it is specified by a PMF $p_X(x_k)=P(X=x_k)$ for a discrete RV, by a PDF $f_X(x)$ with $P(X\in B)=\int_Bf_X(x)\,dx$ for a continuous RV, or in all cases by the CDF $F_X(x)=P(X\le x)$. A **family** of distributions is a formula with adjustable **parameters** ($n,p$ for the binomial; $\mu,\sigma$ for the Gaussian).

### Discrete families

> [!definition] Bernoulli distribution
> $X\sim\operatorname{Bernoulli}(p)$: $P(X=1)=p$, $P(X=0)=1-p$. The indicator of any event with probability $p$; a single bit with error probability $p$.

> [!definition] Binomial distribution
> $X\sim\operatorname{Bin}(n,p)$: the number of successes in $n$ independent Bernoulli($p$) trials,
> $$
> p_X(k)=\binom{n}{k}p^k(1-p)^{n-k},\qquad k=0,1,\dots,n .
> $$

This is the law derived in Section 19.4; Example 19.19 was $\operatorname{Bin}(3,0.1)$. The PMF rises to a peak near $k\approx np$ and then falls; the ratio of successive terms, $p_X(k+1)/p_X(k)=\frac{(n-k)p}{(k+1)(1-p)}$, exceeds 1 exactly when $k<(n+1)p-1$, which locates the mode at $\lfloor(n+1)p\rfloor$.

> [!definition] Geometric distribution
> $N\sim\operatorname{Geom}(p)$: the number of independent trials up to and including the first success,
> $$
> p_N(k)=(1-p)^{k-1}p,\qquad k=1,2,3,\dots
> $$
> with tail $P(N>k)=(1-p)^k$ (the first $k$ trials all fail).

The geometric distribution is the only **memoryless** discrete distribution: $P(N>m+k\mid N>m)=(1-p)^{m+k}/(1-p)^m=P(N>k)$. Having already failed $m$ times does not make success any closer.

> [!definition] Poisson distribution
> $X\sim\operatorname{Poisson}(\lambda)$, $\lambda>0$:
> $$
> p_X(k)=e^{-\lambda}\frac{\lambda^k}{k!},\qquad k=0,1,2,\dots
> $$
> The probabilities sum to $e^{-\lambda}e^{\lambda}=1$ by the exponential series.

> [!theorem] Poisson limit of the binomial
> If $n\to\infty$ and $p\to0$ with $np=\lambda$ fixed, then $\binom{n}{k}p^k(1-p)^{n-k}\to e^{-\lambda}\lambda^k/k!$ for each fixed $k$.

*Proof.* Write $p=\lambda/n$:
$$
\binom{n}{k}\left(\frac{\lambda}{n}\right)^k\left(1-\frac{\lambda}{n}\right)^{n-k}
=\frac{\lambda^k}{k!}\cdot\underbrace{\frac{n(n-1)\cdots(n-k+1)}{n^k}}_{\to1}\cdot\underbrace{\left(1-\frac{\lambda}{n}\right)^{n}}_{\to e^{-\lambda}}\cdot\underbrace{\left(1-\frac{\lambda}{n}\right)^{-k}}_{\to1}. \qquad\square
$$

This is the **law of rare events**: whenever a large number of independent opportunities each produce an event with small probability, the count is approximately Poisson with $\lambda$ = expected number of events. Photon arrivals at a photodiode, electrons crossing a barrier (shot noise), packet arrivals from many independent users, defects on a wafer, and bit errors in a long frame all follow this pattern. If events occur "at random" at an average **rate** $r$ per second, the count in an interval of length $\tau$ is $\operatorname{Poisson}(r\tau)$ — the **Poisson process** of [Ch. 21](21-random-processes.md).

```desmos-graph
left=-0.5; right=10.5; bottom=0; top=0.3
---
(0,0.0282)|blue
(1,0.1211)|blue
(2,0.2335)|blue
(3,0.2668)|blue
(4,0.2001)|blue
(5,0.1029)|blue
(6,0.0368)|blue
(7,0.0090)|blue
(8,0.0014)|blue
(9,0.0001)|blue
(10,0)|blue
(0,0.0498)|red
(1,0.1494)|red
(2,0.2240)|red
(3,0.2240)|red
(4,0.1680)|red
(5,0.1008)|red
(6,0.0504)|red
(7,0.0216)|red
(8,0.0081)|red
(9,0.0027)|red
(10,0.0008)|red
```

*PMFs of $\operatorname{Bin}(10,0.3)$ (blue) and $\operatorname{Poisson}(3)$ (red), both with mean 3. The Poisson PMF is slightly wider and has a tail beyond $k=10$; with only $n=10$ trials and $p=0.3$ the Poisson approximation is rough, but it becomes excellent as $n$ grows with $np$ fixed.*

> [!example] Worked Example 19.22 — Errors in a frame (Poisson approximation)
> A 1000-bit frame is sent over a channel with independent bit errors, $p=0.002$. Find the probability that the frame is error-free and the probability that it has 3 or more errors (more than a double-error-correcting code can handle).
>
> **Step 1 — exact model.** $X\sim\operatorname{Bin}(1000,0.002)$. $P(X=0)=0.998^{1000}=0.13506$.
>
> **Step 2 — Poisson model.** $\lambda=np=2$. $P(X=0)\approx e^{-2}=0.13534$ — within $0.2\,\%$.
>
> **Step 3 — three or more errors.**
> $$
> P(X\ge3)=1-e^{-2}\left(1+2+\frac{2^2}{2!}\right)=1-5e^{-2}=0.32332 .
> $$
> The exact binomial value is $0.32332$ to five decimals. ✓
>
> One frame in three would be uncorrectable; the design must lower $p$, shorten the frame or use a stronger code.

> [!example] Worked Example 19.23 — Retransmissions (geometric)
> An ARQ protocol resends a packet until it is acknowledged; each attempt succeeds independently with probability $0.8$. Find the probability that more than 3 attempts are needed, and the probability that exactly 2 attempts are needed given that the first failed.
>
> **Step 1 — tail.** $P(N>3)=(0.2)^3=0.008$.
>
> **Step 2 — memorylessness.** $P(N=2\mid N>1)=P(N=1)=0.8$ by the memoryless property. Directly: $P(N=2)/P(N>1)=(0.2)(0.8)/0.2=0.8$. ✓

### Continuous families

> [!definition] Uniform distribution
> $X\sim\mathcal{U}(a,b)$: $f_X(x)=\dfrac{1}{b-a}$ for $a<x<b$ and $0$ otherwise. Probability is proportional to length: $P(c<X<d)=(d-c)/(b-a)$ for $a\le c<d\le b$.

The uniform distribution models the phase of an oscillator that is not synchronized to anything ($\Theta\sim\mathcal{U}(0,2\pi)$), and the rounding error of a fine quantizer with step $\Delta$ ($e\sim\mathcal{U}(-\Delta/2,\Delta/2)$, Section 19.9).

> [!definition] Exponential distribution
> $T\sim\operatorname{Exp}(\lambda)$, rate $\lambda>0$:
> $$
> f_T(t)=\lambda e^{-\lambda t}\ (t\ge0),\qquad P(T>t)=e^{-\lambda t}\ (t\ge0).
> $$

> [!theorem] Memoryless property
> $T\sim\operatorname{Exp}(\lambda)$ satisfies $P(T>s+t\mid T>s)=P(T>t)$ for all $s,t\ge0$, and it is the only continuous distribution on $[0,\infty)$ with this property.

*Proof of the property.* $P(T>s+t\mid T>s)=e^{-\lambda(s+t)}/e^{-\lambda s}=e^{-\lambda t}$. (Uniqueness: $G(t)=P(T>t)$ must satisfy $G(s+t)=G(s)G(t)$, whose only monotone solutions are exponentials.) $\square$

A component with exponential lifetime does not age: a unit that has survived 10 000 hours is as good as new. This is realistic for electronic parts during their "useful life" — after infant mortality and before wear-out — where failures are caused by random external stresses rather than accumulated damage. The exponential is also the distribution of the time between events of a Poisson process: if events occur at rate $r$, the waiting time $T$ to the first one satisfies $P(T>t)=P(\text{no events in }[0,t])=e^{-rt}$.

> [!definition] Gaussian (normal) distribution
> $X\sim\mathcal{N}(\mu,\sigma^2)$, with $\mu\in\mathbb{R}$, $\sigma>0$:
> $$
> f_X(x)=\frac{1}{\sqrt{2\pi}\,\sigma}\exp\!\left(-\frac{(x-\mu)^2}{2\sigma^2}\right),\qquad -\infty<x<\infty .
> $$
> The **standard** Gaussian is $Z\sim\mathcal{N}(0,1)$, with PDF $\varphi(z)=\frac{1}{\sqrt{2\pi}}e^{-z^2/2}$, CDF $\Phi(z)$ and tail $Q(z)=1-\Phi(z)$.

*Normalization.* With $I=\int_{-\infty}^{\infty}e^{-z^2/2}dz$, convert $I^2$ to a double integral over the plane and use polar coordinates ([Ch. 6](../part-2-calculus/06-multivariable-calculus.md)):
$$
I^2=\int_{-\infty}^{\infty}\!\int_{-\infty}^{\infty}e^{-(u^2+v^2)/2}\,du\,dv=\int_0^{2\pi}\!\int_0^\infty e^{-r^2/2}\,r\,dr\,d\phi=2\pi ,
$$
so $I=\sqrt{2\pi}$, and the substitution $z=(x-\mu)/\sigma$ handles general $\mu,\sigma$.

The Gaussian is symmetric about $\mu$, has its inflection points at $\mu\pm\sigma$, and its tails fall off faster than any exponential. Its pre-eminence comes from the **central limit theorem** ([Ch. 20](20-statistical-analysis.md)): the sum of many small independent contributions is approximately Gaussian whatever the individual distributions. Thermal noise (the sum of the effects of enormous numbers of electrons), aggregate interference, and manufacturing tolerances that result from many small process variations are all modelled as Gaussian. Gaussian probabilities are computed from tables of $\Phi$ or $Q$ (Section 19.8).

```desmos-graph
left=-6; right=8; bottom=0; top=0.85
---
y=\frac{1}{\sqrt{2\pi}}e^{-\frac{x^{2}}{2}}|blue
y=\frac{1}{0.5\sqrt{2\pi}}e^{-\frac{x^{2}}{2\cdot0.25}}|red
y=\frac{1}{2\sqrt{2\pi}}e^{-\frac{(x-2)^{2}}{8}}|green
```

*Gaussian PDFs: $\mathcal{N}(0,1)$ (blue), $\mathcal{N}(0,0.25)$ (red) and $\mathcal{N}(2,4)$ (green). The mean $\mu$ shifts the curve; the standard deviation $\sigma$ scales its width by $\sigma$ and its height by $1/\sigma$, keeping the area equal to 1.*

> [!definition] Rayleigh distribution
> $R\sim\operatorname{Rayleigh}(\sigma)$:
> $$
> f_R(r)=\frac{r}{\sigma^2}e^{-r^2/(2\sigma^2)}\ (r\ge0),\qquad P(R>r)=e^{-r^2/(2\sigma^2)} .
> $$
> It is the distribution of the envelope $R=\sqrt{X^2+Y^2}$ when $X,Y$ are independent $\mathcal{N}(0,\sigma^2)$ — the amplitude of narrowband Gaussian noise, or of a radio signal arriving over many scattered paths with no direct line of sight. Its mean-square value is $E[R^2]=2\sigma^2$.

> [!example] Worked Example 19.24 — Component lifetime and deep fades
> (a) A power MOSFET has an exponential lifetime with failure rate $\lambda=2\times10^{-5}$ per hour. Find the probability that it survives one year of continuous operation (8760 h), and the probability that it survives a second year given that it survived the first. (b) A Rayleigh-fading signal has mean-square amplitude $E[R^2]=2\sigma^2$. Find the probability that the amplitude is more than 20 dB below its RMS value.
>
> **(a)** $P(T>8760)=e^{-2\times10^{-5}\times8760}=e^{-0.1752}=0.839$. By memorylessness, the conditional probability of surviving the second year is also $0.839$, and of surviving both years $0.839^2=0.704$.
>
> **(b)** 20 dB below RMS in amplitude means $R<0.1\sqrt{2\sigma^2}$, i.e. $R^2<0.01\times2\sigma^2$:
> $$
> P(R<0.1\,R_{\text{rms}})=1-e^{-0.02\sigma^2/(2\sigma^2)}=1-e^{-0.01}=0.00995 .
> $$
> About 1 % of the time the signal is in a 20 dB fade; by the same formula, 10 % of the time it is in a 10 dB fade ($1-e^{-0.1}=0.095$). Each extra 10 dB of fade margin buys a factor of 10 in outage probability — the classic Rayleigh rule of thumb.

```desmos-graph
left=0; right=5; bottom=0; top=2.1
---
y=2e^{-2x}|blue|0<x<5
y=0.5e^{-0.5x}|red|0<x<5
y=xe^{-\frac{x^{2}}{2}}|green|0<x<5
```

*Exponential PDFs with $\lambda=2$ (blue) and $\lambda=0.5$ (red), both maximal at $t=0$: the most likely lifetime is "immediately", though the median is $\ln2/\lambda$. The Rayleigh PDF with $\sigma=1$ (green) vanishes at 0, peaks at $r=\sigma$ and has a Gaussian-like tail.*

### Joint distributions: a first look

Engineering problems usually involve several random quantities at once: the transmitted and received bits, the in-phase and quadrature noise components, the lifetimes of two redundant units. Their *individual* distributions do not say how they interact; that information is in the **joint distribution**. We introduce the essentials for two random variables; correlation and covariance as data-analysis tools are developed in [Ch. 20](20-statistical-analysis.md), and families of RVs indexed by time in [Ch. 21](21-random-processes.md).

> [!definition] Joint PMF, joint PDF, joint CDF
> For two random variables $X,Y$ on the same probability space:
> - the **joint CDF** is $F_{XY}(x,y)=P(X\le x,\ Y\le y)$;
> - if both are discrete, the **joint PMF** is $p_{XY}(x_i,y_k)=P(X=x_i,\ Y=y_k)$, with $\sum_i\sum_kp_{XY}(x_i,y_k)=1$;
> - they are **jointly continuous** with **joint PDF** $f_{XY}(x,y)\ge0$ if, for every region $D$ of the plane,
> $$
> P\bigl((X,Y)\in D\bigr)=\iint_Df_{XY}(x,y)\,dx\,dy,\qquad \iint_{\mathbb{R}^2}f_{XY}\,dx\,dy=1 .
> $$
> Then $f_{XY}(x,y)=\dfrac{\partial^2F_{XY}}{\partial x\,\partial y}$.

From the joint distribution we recover each variable alone — the **marginal** distributions — by summing or integrating out the other variable (law of total probability):

$$
p_X(x_i)=\sum_kp_{XY}(x_i,y_k),\qquad f_X(x)=\int_{-\infty}^{\infty}f_{XY}(x,y)\,dy .
$$

The **conditional** distribution of $Y$ given $X=x$ is the joint divided by the marginal, exactly as in Section 19.2:

$$
p_{Y\mid X}(y_k\mid x_i)=\frac{p_{XY}(x_i,y_k)}{p_X(x_i)},\qquad f_{Y\mid X}(y\mid x)=\frac{f_{XY}(x,y)}{f_X(x)}\quad(f_X(x)>0).
$$

> [!definition] Independent random variables
> $X$ and $Y$ are **independent** if $P(X\in A,\ Y\in B)=P(X\in A)P(Y\in B)$ for all events $A,B$. Equivalently, the joint CDF (or joint PMF, or joint PDF) factors into the product of the marginals:
> $$
> F_{XY}(x,y)=F_X(x)F_Y(y)\quad\text{for all }x,y,\qquad\text{or}\qquad f_{XY}(x,y)=f_X(x)f_Y(y).
> $$

The marginals determine the joint distribution only in the independent case; in general, many different joint distributions share the same marginals.

> [!example] Worked Example 19.25 — Joint PMF of a binary channel
> For the channel of Example 19.11 ($p_0=0.6$, $P(Y=1\mid X=0)=0.05$, $P(Y=0\mid X=1)=0.10$), tabulate the joint PMF of input $X$ and output $Y$, find the marginals, and decide whether $X$ and $Y$ are independent.
>
> **Step 1 — joint PMF** by the multiplication rule $p_{XY}(x,y)=P(Y=y\mid X=x)P(X=x)$:
>
> | | $Y=0$ | $Y=1$ | $p_X$ |
> |---|---|---|---|
> | $X=0$ | $0.57$ | $0.03$ | $0.60$ |
> | $X=1$ | $0.04$ | $0.36$ | $0.40$ |
> | $p_Y$ | $0.61$ | $0.39$ | $1$ |
>
> **Step 2 — marginals** are the row and column sums (written in the margins — hence the name).
>
> **Step 3 — independence test.** $p_X(0)p_Y(0)=0.6\times0.61=0.366\ne0.57$. Not independent — fortunately, since a channel whose output is independent of its input carries no information.

> [!example] Worked Example 19.26 — A joint PDF on a triangle
> $X$ and $Y$ have joint PDF $f_{XY}(x,y)=2$ on the triangle $0<y<x<1$ and $0$ elsewhere. (For example, $X$ and $Y$ could be the larger and smaller of two independent $\mathcal{U}(0,1)$ readings.) Find the marginals, the conditional PDF of $Y$ given $X=x$, and $P(X+Y<1)$.
>
> **Step 1 — check normalization.** The triangle has area $1/2$; $2\times\frac12=1$. ✓
>
> **Step 2 — marginals.**
> $$
> f_X(x)=\int_0^x2\,dy=2x\quad(0<x<1),\qquad f_Y(y)=\int_y^12\,dx=2(1-y)\quad(0<y<1).
> $$
>
> **Step 3 — conditional.** $f_{Y\mid X}(y\mid x)=2/(2x)=1/x$ for $0<y<x$: given $X=x$, $Y$ is uniform on $(0,x)$.
>
> **Step 4 — a probability.** The region $\{y<x,\ x+y<1,\ y>0\}$ is the triangle with vertices $(0,0)$, $(1,0)$, $(\frac12,\frac12)$, of area $\frac14$. So $P(X+Y<1)=2\times\frac14=\frac12$.
>
> $X$ and $Y$ are not independent: $f_X(x)f_Y(y)=4x(1-y)$ is not equal to $2$, and in any case the support (a triangle, not a rectangle) already rules out independence.

### Sums of independent random variables

The distribution of a sum $Z=X+Y$ arises constantly: signal plus noise, total delay of two stages, combined tolerance of two resistors in series, total lifetime of a unit and its cold spare. For independent $X,Y$, the event $\{Z\le z\}$ is the region $x+y\le z$, and

$$
F_Z(z)=\int_{-\infty}^{\infty}\!\int_{-\infty}^{z-x}f_X(x)f_Y(y)\,dy\,dx .
$$

Differentiating with respect to $z$ gives the central result.

> [!theorem] Distribution of a sum of independent RVs
> If $X$ and $Y$ are independent with PDFs $f_X,f_Y$, then $Z=X+Y$ has PDF
> $$
> f_Z(z)=\int_{-\infty}^{\infty}f_X(x)\,f_Y(z-x)\,dx=(f_X*f_Y)(z),
> $$
> the **convolution** of the two densities. For independent integer-valued RVs, $p_Z(n)=\sum_kp_X(k)\,p_Y(n-k)$, the discrete convolution.

This is the same convolution that gives the output of an LTI system ([Ch. 16](../part-5-transforms/16-fourier-transforms.md), [Ch. 17](../part-5-transforms/17-z-transforms.md)), and it can be computed the same way — by transforms. The transform of a density is its characteristic function (Section 19.9), and convolution becomes multiplication.

> [!example] Worked Example 19.27 — Sums of uniform and exponential RVs
> (a) Two independent errors $X,Y\sim\mathcal{U}(0,1)$ add. Find the PDF of $Z=X+Y$. (b) A unit with $\operatorname{Exp}(\lambda)$ lifetime $T_1$ has a cold spare, switched in perfectly on failure, with independent $\operatorname{Exp}(\lambda)$ lifetime $T_2$. Find the reliability of the pair, and evaluate it for $\lambda=2\times10^{-5}\ \text{h}^{-1}$ at one year.
>
> **(a)** $f_X(x)f_Y(z-x)=1$ when $0<x<1$ and $0<z-x<1$, i.e. $\max(0,z-1)<x<\min(1,z)$. The length of that interval is $z$ for $0<z\le1$ and $2-z$ for $1<z<2$:
> $$
> f_Z(z)=\begin{cases}z,&0<z\le1\\2-z,&1<z<2\\0,&\text{otherwise}\end{cases}
> $$
> — a **triangle**, the convolution of two rectangles. Adding more uniforms smooths the shape further towards a Gaussian (central limit theorem).
>
> **(b)** For $z>0$,
> $$
> f_Z(z)=\int_0^z\lambda e^{-\lambda x}\,\lambda e^{-\lambda(z-x)}\,dx=\lambda^2e^{-\lambda z}\int_0^zdx=\lambda^2ze^{-\lambda z},
> $$
> the **Erlang-2** (gamma) density. The reliability is
> $$
> R(t)=P(Z>t)=\int_t^\infty\lambda^2ze^{-\lambda z}dz=(1+\lambda t)e^{-\lambda t}.
> $$
> With $\lambda t=0.1752$: $R=1.1752\,e^{-0.1752}=0.9863$, against $0.8393$ for a single unit and $1-(1-0.8393)^2=0.9742$ for two units running in *active* parallel. A cold spare beats a hot one because it does not age while waiting (assuming a perfect switch).

```desmos-graph
left=-0.5; right=4; bottom=0; top=1.2
---
y=1|blue|0<x<1
y=x|red|0<x<1
y=2-x|red|1<x<2
y=xe^{-x}|green|0<x<4
```

*Convolution in action: the $\mathcal{U}(0,1)$ density (blue) convolved with itself gives the triangle (red); two independent $\operatorname{Exp}(1)$ lifetimes add to the Erlang-2 density $ze^{-z}$ (green), which, unlike the exponential, is zero at the origin — two failures rarely happen in quick succession.*

> [!warning] Common Pitfall — marginals do not determine the joint
> Knowing that $X$ and $Y$ are each $\mathcal{N}(0,1)$ does not determine $P(X+Y>2)$: if $Y=X$ the answer is $Q(1)=0.159$; if $Y=-X$ it is $0$; if they are independent it is $Q(\sqrt2)=0.079$. Any computation involving two random variables needs their joint distribution, and assuming independence is a modelling decision that must be justified.

> [!tip] Engineering Note — independence makes joint models cheap
> Specifying a general joint PMF of $n$ binary variables needs $2^n-1$ numbers; under independence, $n$ numbers suffice. This exponential saving is why independence (or conditional independence given a few hidden variables, as in the fading channel) underlies almost every tractable model of noise, errors and failures.

## 19.7 Probability density and mass functions

### Motivation

The PMF and PDF are the working descriptions of a random variable: they are what you sketch, fit to histograms, multiply by likelihoods and integrate against. The PMF is simple — a list of probabilities. The PDF needs more care, because its values are *not* probabilities, it can exceed 1, it can contain impulses, and it transforms in a non-obvious way when the variable is passed through a nonlinearity. This section states the properties precisely, shows how to treat discrete, continuous and mixed variables in one notation, and derives the transformation rule that gives the density at the output of a memoryless system.

### The probability mass function

> [!definition] PMF
> The **probability mass function** of a discrete random variable $X$ with range $S_X=\{x_1,x_2,\dots\}$ is
> $$
> p_X(x)=P(X=x).
> $$
> It satisfies $p_X(x)\ge0$, $p_X(x)=0$ for $x\notin S_X$, and $\displaystyle\sum_{k}p_X(x_k)=1$. For any set $B$, $P(X\in B)=\sum_{x_k\in B}p_X(x_k)$.

Any non-negative sequence summing to 1 is a valid PMF, and constructing a model often means choosing such a sequence from a family and fixing its parameters.

### The probability density function

> [!definition] PDF
> A random variable $X$ is **continuous** with **probability density function** $f_X$ if $f_X(x)\ge0$, $\int_{-\infty}^{\infty}f_X(x)\,dx=1$, and for all $a<b$
> $$
> P(a<X\le b)=\int_a^bf_X(x)\,dx .
> $$

The name "density" is literal. For a small interval of width $\Delta x$ around a point where $f_X$ is continuous,

$$
P(x<X\le x+\Delta x)\approx f_X(x)\,\Delta x ,
$$

so $f_X(x)$ is probability *per unit of $x$*. If $X$ is a voltage, $f_X$ has units of $\text{V}^{-1}$. This has three consequences that trip up newcomers:

1. **$f_X(x)$ is not a probability.** It can exceed 1: the $\mathcal{U}(0,0.1)$ density equals 10 on its support.
2. **Only areas are probabilities.** $P(X=x)=\int_x^xf_X=0$ for every $x$.
3. **Densities change under a change of units.** If $Y=1000X$ (the same voltage in mV), then $f_Y(y)=f_X(y/1000)/1000$, so that $f_Y(y)\,dy=f_X(x)\,dx$ (see the transformation rule below).

A histogram of measured values, scaled so that its total area is 1 (bar height $=$ fraction of samples in the bin divided by bin width), approximates the PDF; that is how densities are estimated from data in [Ch. 20](20-statistical-analysis.md).

Useful single-number descriptions read off a PDF or PMF are the **mode** (the location of the maximum — the "most likely" value), the **median** $m$ (with $P(X\le m)=1/2$, Section 19.8), and the **mean** (the balance point, Section 19.9). For a symmetric unimodal density like the Gaussian, all three coincide; for a skewed one like the exponential they differ (mode $0$, median $\ln2/\lambda$, mean $1/\lambda$).

> [!example] Worked Example 19.28 — Normalizing a density
> The normalized position $X$ of a laser spot on a detector has PDF $f_X(x)=c\,x(1-x)$ for $0<x<1$ and $0$ otherwise. Find $c$, $P(X>0.5)$ and $P(0.25<X<0.75)$.
>
> **Step 1 — normalize.** $\int_0^1x(1-x)\,dx=\frac12-\frac13=\frac16$, so $c=6$.
>
> **Step 2 — by symmetry** about $x=\frac12$, $P(X>0.5)=0.5$.
>
> **Step 3 — integrate.** With antiderivative $G(x)=3x^2-2x^3$,
> $$
> P(0.25<X<0.75)=G(0.75)-G(0.25)=(1.6875-0.84375)-(0.1875-0.03125)=0.84375-0.15625=0.6875 .
> $$
> The maximum value of the density is $f_X(0.5)=1.5$ — larger than 1, which is perfectly legitimate.

### Discrete and mixed variables as densities: Dirac impulses

Using the Dirac delta $\delta(x)$ ([Ch. 16](../part-5-transforms/16-fourier-transforms.md)), a discrete random variable can be given a *generalized* density

$$
f_X(x)=\sum_kp_X(x_k)\,\delta(x-x_k),
$$

because $\int_a^bf_X(x)\,dx$ then picks up exactly the masses $p_X(x_k)$ with $a<x_k\le b$. A mixed random variable has a density that is the sum of an ordinary part and impulses. The limiter output of Example 19.21, for instance, has

$$
f_Y(y)=0.25\,\delta(y+1)+0.25\,\delta(y-1)+0.25\,\operatorname{rect}(y/2),
$$

where $\operatorname{rect}(y/2)=1$ for $\lvert y\rvert<1$. This unified notation lets every formula of this chapter — expectations, transformations, convolutions — be written once for all three types. The CDF (Section 19.8) achieves the same unification without generalized functions.

### Conditional densities given an event

If $B$ is an event with $P(B)>0$, the **conditional PDF** of $X$ given $B$ is the derivative of the conditional CDF $P(X\le x\mid B)$. When $B=\{X\in D\}$ is itself defined by $X$, the result is the original density restricted to $D$ and renormalized:

$$
f_{X\mid B}(x)=\begin{cases}\dfrac{f_X(x)}{P(X\in D)},&x\in D\\[2mm]0,&x\notin D.\end{cases}
$$

For an exponential lifetime and $B=\{T>s\}$ this gives $f_{T\mid B}(t)=\lambda e^{-\lambda(t-s)}$ for $t>s$: the same exponential shape shifted to start at $s$ — memorylessness again.

### Bayes' theorem with densities

When the observation is continuous, the event $\{X=x\}$ has probability zero and Section 19.3's formula cannot be applied literally. Conditioning on the small interval $\{x<X\le x+\Delta x\}$ and letting $\Delta x\to0$, the factors $\Delta x$ cancel, leaving

$$
P(H_i\mid X=x)=\frac{f_{X\mid H_i}(x)\,P(H_i)}{\sum_jf_{X\mid H_j}(x)\,P(H_j)} ,
$$

i.e. Bayes' theorem with *conditional densities* playing the role of likelihoods. The MAP rule compares $f_{X\mid H_i}(x)P(H_i)$.

> [!example] Worked Example 19.29 — MAP threshold for a binary signal in Gaussian noise
> A receiver observes $X=s+N$, where $s=-A$ (symbol 0) with probability $p_0$ or $s=+A$ (symbol 1) with probability $p_1$, and $N\sim\mathcal{N}(0,\sigma^2)$ is independent noise. Find the MAP decision threshold, and evaluate it and the error probability for $A=1$ V, $\sigma=0.5$ V, $p_0=0.75$, $p_1=0.25$.
>
> **Step 1 — likelihoods.** Given symbol 1, $X\sim\mathcal{N}(A,\sigma^2)$; given symbol 0, $X\sim\mathcal{N}(-A,\sigma^2)$.
>
> **Step 2 — likelihood-ratio test.** Decide 1 when $p_1f_{X\mid1}(x)>p_0f_{X\mid0}(x)$:
> $$
> \frac{f_{X\mid1}(x)}{f_{X\mid0}(x)}=\exp\!\left(\frac{(x+A)^2-(x-A)^2}{2\sigma^2}\right)=\exp\!\left(\frac{2Ax}{\sigma^2}\right)>\frac{p_0}{p_1}.
> $$
> Taking logarithms: decide 1 when $x>x_{\text{th}}$, with
> $$
> x_{\text{th}}=\frac{\sigma^2}{2A}\ln\frac{p_0}{p_1}.
> $$
>
> **Step 3 — numbers.** $x_{\text{th}}=\dfrac{0.25}{2}\ln3=0.1373$ V. The threshold moves *towards the less likely symbol*, so that more ambiguous observations are attributed to the more likely one.
>
> **Step 4 — error probability.** An error occurs if symbol 0 is sent and $X>x_{\text{th}}$, or symbol 1 is sent and $X<x_{\text{th}}$:
> $$
> P_e=p_0\,Q\!\left(\frac{x_{\text{th}}+A}{\sigma}\right)+p_1\,Q\!\left(\frac{A-x_{\text{th}}}{\sigma}\right)=0.75\,Q(2.275)+0.25\,Q(1.725)=0.01916 .
> $$
> With the "natural" threshold $0$ (the ML rule), $P_e=Q(A/\sigma)=Q(2)=0.02275$. The MAP threshold reduces the error probability by about $16\,\%$. With equal priors, $\ln(p_0/p_1)=0$ and the two rules coincide.

```desmos-graph
left=-3; right=3; bottom=0; top=0.65
---
y=0.75\cdot\frac{1}{0.5\sqrt{2\pi}}e^{-\frac{(x+1)^{2}}{0.5}}|blue
y=0.25\cdot\frac{1}{0.5\sqrt{2\pi}}e^{-\frac{(x-1)^{2}}{0.5}}|red
x=0.1373|black|dashed
```

*Prior-weighted likelihoods $p_0f_{X\mid0}(x)$ (blue) and $p_1f_{X\mid1}(x)$ (red) for Example 19.29. The MAP threshold (dashed) is where the two curves cross, at $x=0.137$ V rather than midway at 0; the error probability is the area under the lower curve on each side of the threshold.*

### Functions of a random variable: the density transformation

We now find the PDF of $Y=g(X)$ when $X$ is continuous. Start with the simplest case.

> [!theorem] Transformation by a monotone function
> Let $X$ have PDF $f_X$ and let $g$ be differentiable and strictly monotone (increasing or decreasing) on the range of $X$, with inverse $x=g^{-1}(y)$. Then $Y=g(X)$ has PDF
> $$
> f_Y(y)=f_X\bigl(g^{-1}(y)\bigr)\left\lvert\frac{dx}{dy}\right\rvert=\frac{f_X(x)}{\lvert g'(x)\rvert}\bigg|_{x=g^{-1}(y)} .
> $$

*Derivation.* For increasing $g$, $\{Y\le y\}=\{X\le g^{-1}(y)\}$, so $F_Y(y)=F_X(g^{-1}(y))$; differentiate by the chain rule to get $f_X(g^{-1}(y))\,\frac{d}{dy}g^{-1}(y)$. For decreasing $g$, $\{Y\le y\}=\{X\ge g^{-1}(y)\}$, so $F_Y(y)=1-F_X(g^{-1}(y))$, whose derivative carries a minus sign that the absolute value absorbs. Physically: probability is conserved, $f_Y(y)\lvert dy\rvert=f_X(x)\lvert dx\rvert$ — a steep part of $g$ stretches a small $dx$ into a large $dy$ and thins the density. $\square$

> [!theorem] Transformation by a general function
> If the equation $y=g(x)$ has finitely many real roots $x_1(y),\dots,x_m(y)$, at each of which $g'(x_i)\ne0$, then
> $$
> f_Y(y)=\sum_{i=1}^{m}\frac{f_X(x_i)}{\lvert g'(x_i)\rvert} ,
> $$
> and $f_Y(y)=0$ if there are no roots. Intervals on which $g$ is constant produce point masses (impulses) in $f_Y$, as in the limiter.

*Reason.* The event $\{y<Y\le y+dy\}$ is the union of disjoint small intervals of $x$, one around each root, of widths $dy/\lvert g'(x_i)\rvert$; add their probabilities.

> [!example] Worked Example 19.30 — Linear transformation of a Gaussian
> A sensor output $X\sim\mathcal{N}(\mu,\sigma^2)$ is amplified and offset: $Y=aX+b$ with $a\ne0$. Find the distribution of $Y$, and the density at the output of a gain-10, offset $+0.5$ V stage when $X\sim\mathcal{N}(0.1,\,0.02^2)$ V.
>
> **Step 1 — inverse.** $x=(y-b)/a$, $\lvert dx/dy\rvert=1/\lvert a\rvert$.
>
> **Step 2 — substitute.**
> $$
> f_Y(y)=\frac{1}{\lvert a\rvert}\cdot\frac{1}{\sqrt{2\pi}\sigma}\exp\!\left(-\frac{\bigl(\frac{y-b}{a}-\mu\bigr)^2}{2\sigma^2}\right)=\frac{1}{\sqrt{2\pi}\,\lvert a\rvert\sigma}\exp\!\left(-\frac{(y-a\mu-b)^2}{2a^2\sigma^2}\right).
> $$
> So $Y\sim\mathcal{N}(a\mu+b,\ a^2\sigma^2)$: **a linear function of a Gaussian is Gaussian.**
>
> **Step 3 — numbers.** $Y\sim\mathcal{N}(10(0.1)+0.5,\ (10\times0.02)^2)=\mathcal{N}(1.5,\ 0.2^2)$ V. In particular $Z=(X-\mu)/\sigma$ is $\mathcal{N}(0,1)$ — **standardization**, which reduces every Gaussian probability to the standard tables.

> [!example] Worked Example 19.31 — Square-law detector
> Gaussian noise $X\sim\mathcal{N}(0,\sigma^2)$ drives a square-law detector, $Y=X^2$. Find $f_Y$ and the probability that the detector output exceeds $4\sigma^2$.
>
> **Step 1 — roots.** For $y>0$, $x^2=y$ has roots $x_{1,2}=\pm\sqrt y$, with $\lvert g'(x)\rvert=2\lvert x\rvert=2\sqrt y$. For $y<0$ there are no roots.
>
> **Step 2 — add contributions.**
> $$
> f_Y(y)=\frac{f_X(\sqrt y)+f_X(-\sqrt y)}{2\sqrt y}=\frac{1}{\sqrt{2\pi y}\,\sigma}\,e^{-y/(2\sigma^2)},\qquad y>0 .
> $$
> This is the **chi-square distribution with one degree of freedom** (scaled by $\sigma^2$). It is infinite at $y=0$ — yet integrable, since $y^{-1/2}$ has a finite integral near 0.
>
> **Step 3 — tail probability.** Easier through $X$ than through $f_Y$: $P(Y>4\sigma^2)=P(\lvert X\rvert>2\sigma)=2Q(2)=0.0455$.

> [!example] Worked Example 19.32 — A sinusoid with random phase
> An oscillator's output sampled at a random instant is $Y=A\cos\Theta$, with $\Theta\sim\mathcal{U}(0,2\pi)$. Find the PDF of $Y$ and the fraction of time the sinusoid spends within half its amplitude of zero.
>
> **Step 1 — roots.** For $\lvert y\rvert<A$, $A\cos\theta=y$ has two roots in $(0,2\pi)$, $\theta_1=\arccos(y/A)$ and $\theta_2=2\pi-\theta_1$. At each, $\lvert g'(\theta)\rvert=A\lvert\sin\theta\rvert=\sqrt{A^2-y^2}$.
>
> **Step 2 — transform.** $f_\Theta=1/(2\pi)$ on $(0,2\pi)$, so
> $$
> f_Y(y)=2\cdot\frac{1/(2\pi)}{\sqrt{A^2-y^2}}=\frac{1}{\pi\sqrt{A^2-y^2}},\qquad \lvert y\rvert<A .
> $$
> This is the **arcsine** density. It is smallest at $y=0$ and infinite (but integrable) at $y=\pm A$: a sinusoid moves fastest through zero and lingers near its peaks.
>
> **Step 3 — probability.** $P(\lvert Y\rvert<A/2)=\int_{-A/2}^{A/2}\frac{dy}{\pi\sqrt{A^2-y^2}}=\frac{2}{\pi}\arcsin\frac12=\frac13$.
>
> Only a third of the time is spent in the middle half of the range. This is why a histogram of ADC codes for a full-scale sine input (the "code density test") is bathtub-shaped, and why departures from the arcsine shape reveal the converter's non-linearity.

```desmos-graph
left=-1.3; right=1.3; bottom=0; top=2
---
y=\frac{1}{\pi\sqrt{1-x^{2}}}|blue|-1<x<1
y=0.5|red|-1<x<1
```

*Arcsine density of $\cos\Theta$ with uniform phase (blue), compared with the uniform density on $(-1,1)$ (red). The sinusoid spends most of its time near its extremes; a triangle wave of the same amplitude would give the flat red density.*

> [!warning] Common Pitfall — forgetting the Jacobian
> The most common error in transforming densities is to substitute $x=g^{-1}(y)$ into $f_X$ and stop, omitting the factor $\lvert dx/dy\rvert$. The result is not even normalized. A quick check: integrate your $f_Y$ — it must give 1. A second check: the units of $f_Y$ must be the reciprocal of the units of $Y$.

> [!tip] Engineering Note — memoryless nonlinearities
> The transformation rule is the probabilistic counterpart of a transfer characteristic. Diodes, square-law and envelope detectors, logarithmic amplifiers, companders, limiters and quantizers all map an input distribution into an output distribution by exactly this rule. Designers use it in reverse, too: a compander is chosen so that the output of a speech signal with a peaked (Laplacian-like) amplitude distribution is spread more evenly across the quantizer levels.

## 19.8 Cumulative distributions

### Motivation

The PMF works only for discrete variables and the PDF only for continuous ones (unless we admit impulses). The **cumulative distribution function** works for every random variable, discrete, continuous or mixed, without any generalized functions. It is also the function that answers the questions engineers most often ask — "what fraction of resistors are below 990 Ω?", "what is the probability that the component has failed by time $t$?", "what threshold is exceeded only 0.1 % of the time?" — and it is the natural tool for deriving distributions of functions of random variables, extremes of several variables, and for generating random numbers in simulation.

### Definition and properties

> [!definition] Cumulative distribution function (CDF)
> The **CDF** of a random variable $X$ is
> $$
> F_X(x)=P(X\le x),\qquad -\infty<x<\infty .
> $$

> [!theorem] Properties of the CDF
> Every CDF satisfies:
> 1. $0\le F_X(x)\le1$;
> 2. $F_X$ is **non-decreasing**: $a<b\Rightarrow F_X(a)\le F_X(b)$;
> 3. $\lim_{x\to-\infty}F_X(x)=0$ and $\lim_{x\to+\infty}F_X(x)=1$;
> 4. $F_X$ is **right-continuous**: $\lim_{h\downarrow0}F_X(x+h)=F_X(x)$;
> 5. $P(a<X\le b)=F_X(b)-F_X(a)$;
> 6. $P(X=a)=F_X(a)-F_X(a^-)$, the size of the jump at $a$, where $F_X(a^-)=\lim_{h\downarrow0}F_X(a-h)$.
>
> Conversely, any function with properties 2–4 is the CDF of some random variable.

*Proof sketch.* (2) $\{X\le a\}\subseteq\{X\le b\}$, then monotonicity of $P$. (5) $\{X\le b\}$ is the disjoint union of $\{X\le a\}$ and $\{a<X\le b\}$. (3) and (4) follow from continuity of probability (Section 19.1): the events $\{X\le -n\}$ decrease to $\varnothing$, the events $\{X\le n\}$ increase to $\Omega$, and $\{X\le x+1/n\}$ decrease to $\{X\le x\}$. (6) The events $\{a-1/n<X\le a\}$ decrease to $\{X=a\}$, and their probabilities are $F_X(a)-F_X(a-1/n)$. $\square$

Property 6 is the key to reading a CDF. Wherever $F_X$ is continuous, $P(X=x)=0$; wherever it jumps, the jump is a point mass. Hence:

| Type of RV | Shape of the CDF | Link to PMF/PDF |
|---|---|---|
| Discrete | staircase, flat between jumps | $F_X(x)=\sum_{x_k\le x}p_X(x_k)$; jumps $=p_X(x_k)$ |
| Continuous | continuous (no jumps) | $F_X(x)=\int_{-\infty}^xf_X(u)\,du$, $f_X(x)=F_X'(x)$ |
| Mixed | ramps plus jumps | jumps $=$ point masses; slope $=$ continuous density |

For a continuous RV the relation $f_X=F_X'$ holds wherever $F_X$ is differentiable; for a discrete or mixed RV it holds in the generalized sense, the derivative of a jump of height $p$ being $p\,\delta(x-x_k)$.

> [!warning] Common Pitfall — $<$ versus $\le$
> With the convention $F_X(x)=P(X\le x)$, the probability of an interval depends on whether its endpoints are included when $X$ has point masses there: $P(a\le X\le b)=F_X(b)-F_X(a^-)$, $P(a<X<b)=F_X(b^-)-F_X(a)$. For a discrete count, $P(X<3)=F_X(2)$, not $F_X(3)$. For continuous variables the distinction disappears.

> [!example] Worked Example 19.33 — Staircase CDF of an error count
> For the number of errors $X\sim\operatorname{Bin}(3,0.1)$ of Example 19.19, write $F_X(x)$, and use it to find $P(1\le X\le2)$ and $P(X>1)$.
>
> **Step 1 — accumulate the PMF** $(0.729,0.243,0.027,0.001)$:
> $$
> F_X(x)=\begin{cases}0,&x<0\\0.729,&0\le x<1\\0.972,&1\le x<2\\0.999,&2\le x<3\\1,&x\ge3.\end{cases}
> $$
>
> **Step 2 — interval with point masses.** $P(1\le X\le2)=F_X(2)-F_X(1^-)=0.999-0.729=0.270$. Check: $0.243+0.027=0.270$. ✓
>
> **Step 3 — tail.** $P(X>1)=1-F_X(1)=1-0.972=0.028$.

```desmos-graph
left=-1; right=4; bottom=-0.05; top=1.1
---
y=0|blue|-1<x<0
y=0.729|blue|0<x<1
y=0.972|blue|1<x<2
y=0.999|blue|2<x<3
y=1|blue|3<x<4
(0,0.729)|blue
(1,0.972)|blue
(2,0.999)|blue
(3,1)|blue
```

*Staircase CDF of $\operatorname{Bin}(3,0.1)$. Each jump equals the PMF value at that point; the filled dots mark the right-continuity convention — the CDF takes the upper value at each jump.*

### The Gaussian CDF and the $Q$-function

The Gaussian CDF has no closed form; it is tabulated, or computed from the **error function** $\operatorname{erf}(x)=\frac{2}{\sqrt\pi}\int_0^xe^{-t^2}dt$, available in every numerical library:

$$
\Phi(z)=\int_{-\infty}^{z}\frac{1}{\sqrt{2\pi}}e^{-u^2/2}\,du=\frac12\left[1+\operatorname{erf}\!\left(\frac{z}{\sqrt2}\right)\right],\qquad Q(z)=1-\Phi(z)=\frac12\operatorname{erfc}\!\left(\frac{z}{\sqrt2}\right).
$$

By symmetry $\Phi(-z)=1-\Phi(z)$ and $Q(-z)=1-Q(z)$. For $X\sim\mathcal{N}(\mu,\sigma^2)$, standardization (Example 19.30) gives

$$
F_X(x)=\Phi\!\left(\frac{x-\mu}{\sigma}\right),\qquad P(X>x)=Q\!\left(\frac{x-\mu}{\sigma}\right).
$$

Communications engineers use $Q$ because error probabilities are tail probabilities; statisticians tabulate $\Phi$. A short table covers most needs:

| $z$ | $\Phi(z)$ | $Q(z)$ |
|---|---|---|
| $0$ | $0.5000$ | $0.5000$ |
| $1$ | $0.8413$ | $0.1587$ |
| $1.645$ | $0.9500$ | $0.0500$ |
| $1.960$ | $0.9750$ | $0.0250$ |
| $2$ | $0.9772$ | $0.02275$ |
| $2.326$ | $0.9900$ | $0.0100$ |
| $2.576$ | $0.9950$ | $0.0050$ |
| $3$ | $0.99865$ | $1.350\times10^{-3}$ |
| $4$ | — | $3.167\times10^{-5}$ |
| $5$ | — | $2.867\times10^{-7}$ |
| $6$ | — | $9.866\times10^{-10}$ |

The probabilities of falling within $1$, $2$, $3$ standard deviations of the mean are $0.6827$, $0.9545$ and $0.9973$ — the "68–95–99.7 rule". For large $z$ the bound

$$
\frac{z}{1+z^2}\,\frac{e^{-z^2/2}}{\sqrt{2\pi}}<Q(z)<\frac{1}{z}\,\frac{e^{-z^2/2}}{\sqrt{2\pi}}\qquad(z>0)
$$

shows that $Q$ falls off essentially like $e^{-z^2/2}$; the simpler Chernoff bound $Q(z)\le\frac12e^{-z^2/2}$ is often enough for design estimates.

```desmos-graph
left=-4; right=4; bottom=-0.05; top=1.1
---
y=\frac{1}{2}\left(1+\operatorname{erf}\left(\frac{x}{\sqrt{2}}\right)\right)|blue
y=\frac{1}{2}\left(1-\operatorname{erf}\left(\frac{x}{\sqrt{2}}\right)\right)|red|dashed
y=\frac{1}{\sqrt{2\pi}}e^{-\frac{x^{2}}{2}}|grey|dotted
```

*Standard Gaussian CDF $\Phi(z)$ (blue), $Q$-function $1-\Phi(z)$ (red, dashed) and the PDF $\varphi(z)$ (grey, dotted). The CDF's slope is the PDF: steepest at $z=0$, where the density is largest.*

> [!example] Worked Example 19.34 — Resistor tolerance
> A batch of nominal $1\ \text{k}\Omega$ resistors has values $R\sim\mathcal{N}(1000,\,10^2)\ \Omega$. (a) What fraction lies between $990$ and $1020\ \Omega$? (b) What fraction lies outside a $\pm2\,\%$ tolerance band? (c) What symmetric band contains $99.9\,\%$ of the batch?
>
> **(a)** Standardize: $z_1=(990-1000)/10=-1$, $z_2=(1020-1000)/10=2$.
> $$
> P(990<R<1020)=\Phi(2)-\Phi(-1)=0.97725-0.15866=0.8186 .
> $$
>
> **(b)** $\pm2\,\%$ is $\pm20\ \Omega=\pm2\sigma$: $P(\lvert R-1000\rvert>20)=2Q(2)=0.0455$, i.e. $4.55\,\%$ of parts fail a $2\,\%$ specification.
>
> **(c)** We need $2Q(z)=0.001$, $Q(z)=0.0005$, $z=3.291$: the band is $1000\pm32.9\ \Omega$, i.e. about $\pm3.3\,\%$.

### Quantiles

> [!definition] Quantile, median, percentile
> For $0<u<1$, the **$u$-quantile** of $X$ is
> $$
> x_u=F_X^{-1}(u)=\inf\{x: F_X(x)\ge u\},
> $$
> the smallest $x$ at which the CDF reaches $u$. The **median** is $x_{0.5}$; the **quartiles** are $x_{0.25},x_{0.75}$; the $100u$-th **percentile** is $x_u$.

For a continuous, strictly increasing CDF, $F_X^{-1}$ is the ordinary inverse. The "inf" in the definition handles flat parts and jumps. Quantiles are how specifications are written: "the 99th-percentile latency is 4 ms", "the $B_{10}$ life of a bearing (10 % failed)", "the 1 % outage fade depth". For the exponential distribution, solving $1-e^{-\lambda t}=u$ gives $t_u=-\ln(1-u)/\lambda$; the median life is $\ln2/\lambda$, about $69\,\%$ of the mean life $1/\lambda$.

### Reliability function and hazard rate

For a lifetime $T\ge0$, the complement of the CDF is the **reliability** (or **survival**) function, and its relative rate of decrease is the **hazard rate**:

$$
R(t)=P(T>t)=1-F_T(t),\qquad h(t)=\frac{f_T(t)}{R(t)}=-\frac{d}{dt}\ln R(t).
$$

The hazard rate is the conditional failure density of a unit that has survived to $t$: $P(t<T\le t+dt\mid T>t)\approx h(t)\,dt$. Integrating, $R(t)=\exp\!\left(-\int_0^th(u)\,du\right)$. A constant hazard $h(t)=\lambda$ gives the exponential distribution — the mathematical statement of memorylessness. An increasing hazard means wear-out (electrolytic capacitors, fans, electromigration); a decreasing hazard means infant mortality (latent manufacturing defects, removed by burn-in). The sequence decreasing–constant–increasing over a product's life is the well-known **bathtub curve**; the Weibull family that models each phase is covered in [Ch. 20](20-statistical-analysis.md).

### The CDF method for functions of random variables

The most reliable way to find the distribution of $Y=g(X)$ — always valid, including for non-monotone $g$, flat regions, and mixed outputs — is to compute its CDF directly from the definition:

$$
F_Y(y)=P(g(X)\le y)=P\bigl(X\in\{x: g(x)\le y\}\bigr),
$$

and then differentiate if a density is wanted. The density-transformation formulas of Section 19.7 are shortcuts derived this way.

> [!example] Worked Example 19.35 — Half-wave rectified noise (a mixed CDF)
> Noise $X\sim\mathcal{N}(0,1)$ V passes through an ideal half-wave rectifier, $Y=\max(X,0)$. Find $F_Y$ and describe the distribution.
>
> **Step 1 — $y<0$.** $Y\ge0$ always, so $F_Y(y)=0$.
>
> **Step 2 — $y\ge0$.** $\{Y\le y\}=\{X\le y\}$ (if $X\le0$ then $Y=0\le y$; if $X>0$ then $Y=X$). So $F_Y(y)=\Phi(y)$.
>
> **Step 3 — interpret.** $F_Y$ jumps from $0$ to $\Phi(0)=\frac12$ at $y=0$: a point mass $P(Y=0)=\frac12$ (all negative inputs map to 0). For $y>0$ the CDF follows $\Phi$, with density $\varphi(y)$. In generalized-density form, $f_Y(y)=\frac12\delta(y)+\varphi(y)u(y)$.

```desmos-graph
left=-3; right=3; bottom=-0.05; top=1.1
---
y=0|blue|-3<x<0
y=\frac{1}{2}\left(1+\operatorname{erf}\left(\frac{x}{\sqrt{2}}\right)\right)|blue|0<x<3
(0,0.5)|blue
(0,0)|blue|label:jump of 1/2 at y = 0
```

*CDF of half-wave rectified Gaussian noise: zero for $y<0$, a jump of $1/2$ at $y=0$ (the point mass of all clipped negative inputs), then the Gaussian CDF. A mixed random variable needs no impulses in this description.*

### Extremes: the minimum and the maximum

The CDF method is particularly efficient for the largest and smallest of several independent random variables, which appear as the lifetimes of parallel and series systems and as peak values of noise.

> [!theorem] Distribution of the maximum and minimum
> Let $X_1,\dots,X_n$ be independent with CDFs $F_1,\dots,F_n$. Then
> $$
> F_{\max}(x)=P(\max_iX_i\le x)=\prod_{i=1}^nF_i(x),\qquad P(\min_iX_i>x)=\prod_{i=1}^n\bigl(1-F_i(x)\bigr).
> $$

*Proof.* The maximum is at most $x$ exactly when every $X_i\le x$; the minimum exceeds $x$ exactly when every $X_i>x$. Independence turns the intersections into products. $\square$

For lifetimes these are the reliability formulas of Section 19.4 in CDF language: a series system fails at the *first* component failure, $T_{\text{sys}}=\min_iT_i$, so $R_{\text{sys}}(t)=\prod R_i(t)$; a parallel system fails at the *last*, $T_{\text{sys}}=\max_iT_i$, so $F_{\text{sys}}(t)=\prod F_i(t)$. With exponential components, $\min_iT_i\sim\operatorname{Exp}(\lambda_1+\cdots+\lambda_n)$.

> [!example] Worked Example 19.36 — Peaks of noise
> A comparator sees 1000 independent samples of Gaussian noise $X_i\sim\mathcal{N}(0,\sigma^2)$ per frame. Its threshold is set at $3\sigma$, "because noise exceeds $3\sigma$ only $0.13\,\%$ of the time". Find the probability of at least one false trigger per frame, and the threshold needed to make it $1\,\%$.
>
> **Step 1 — maximum.** $P(\max_iX_i\le3\sigma)=\Phi(3)^{1000}=0.998650^{1000}=0.2590$.
>
> **Step 2 — false trigger.** $P(\text{at least one})=1-0.2590=0.741$. Almost three frames in four trigger falsely!
>
> **Step 3 — redesign.** We need $\Phi(z)^{1000}=0.99$, i.e. $\Phi(z)=0.99^{1/1000}=0.99998995$, so $Q(z)=1.005\times10^{-5}$ and $z\approx4.26$. The threshold must be about $4.3\sigma$.
>
> When many independent chances are taken, a rare event per trial becomes a likely event per frame (compare the union bound $1000\times0.00135=1.35$, which is useless here because it exceeds 1). This is the "look-elsewhere" effect that governs radar CFAR thresholds and test-limit settings on automatic test equipment.

### Generating random variables: the inverse-transform method

The quantile function also runs the other way. If $U\sim\mathcal{U}(0,1)$ and $F$ is any CDF, then

$$
X=F^{-1}(U)\quad\text{has CDF } F,
$$

because $P(F^{-1}(U)\le x)=P(U\le F(x))=F(x)$. Conversely, if $X$ is continuous then $F_X(X)\sim\mathcal{U}(0,1)$ (the **probability integral transform**). This is how simulations turn a uniform random-number generator into samples from any distribution: for $\operatorname{Exp}(\lambda)$, $X=-\ln(1-U)/\lambda$; for Rayleigh($\sigma$), $X=\sigma\sqrt{-2\ln(1-U)}$; for a discrete distribution, choose the smallest $x_k$ with $F(x_k)\ge U$. Monte Carlo methods built on this idea are developed in [Ch. 24](../part-7-numerical-computational/24-computational-mathematics.md).

> [!example] Worked Example 19.37 — Simulating component lifetimes
> A reliability simulation needs lifetimes with $\operatorname{Exp}(\lambda)$, $\lambda=2\times10^{-5}\ \text{h}^{-1}$. The uniform generator returns $u=0.5$ and $u=0.95$. What lifetimes are produced, and what do they mean?
>
> **Step 1 — formula.** $t=-\ln(1-u)/\lambda$.
>
> **Step 2 — numbers.** $u=0.5$: $t=\ln2/(2\times10^{-5})=34\,657$ h, the median life. $u=0.95$: $t=-\ln(0.05)/(2\times10^{-5})=149\,787$ h, the 95th-percentile life.
>
> **Step 3 — check.** $F(34\,657)=1-e^{-0.6931}=0.5$. ✓ Each uniform sample is mapped to the lifetime with that cumulative probability.

> [!tip] Engineering Note — read specifications as CDF points
> Nearly every probabilistic specification is a point on a CDF or its complement: a bit-error rate is $P(\text{error})$, an outage probability is $F_{\text{SNR}}(\text{threshold})$, a "six-sigma" process puts the specification limit where $Q(z)$ is tiny, an MTTF is a mean but a warranty period is a quantile. Plotting the empirical CDF of measured data (Ch. 20) is the quickest way to compare a batch against such a specification.

## 19.9 Expectation and variance

### Motivation

A full distribution is often more information than a design needs. An engineer wants to know the DC level of a noisy signal, its RMS fluctuation, the average number of retransmissions, the mean time to failure, the noise power at the output of an ADC. These are all **expectations** — probability-weighted averages — of the random variable or of a function of it. The mean locates the distribution, the variance measures its spread, higher moments describe its shape, and the moment-generating and characteristic functions package all the moments into a single transform that turns the convolution of densities into multiplication. This section develops these tools and the inequalities that make them useful even when the distribution itself is unknown.

### Expected value

> [!definition] Expectation (mean)
> The **expected value** (mean) of $X$ is
> $$
> E[X]=\mu_X=\begin{cases}\displaystyle\sum_kx_k\,p_X(x_k),&X\text{ discrete},\\[3mm]\displaystyle\int_{-\infty}^{\infty}x\,f_X(x)\,dx,&X\text{ continuous},\end{cases}
> $$
> provided the sum or integral converges absolutely ($\sum\lvert x_k\rvert p_X(x_k)<\infty$ or $\int\lvert x\rvert f_X\,dx<\infty$). Otherwise the mean does not exist.

The mean is the **centre of mass** of the probability distribution: if unit mass is distributed along the $x$-axis according to $f_X$, $E[X]$ is the balance point. Its link with experiment is the **law of large numbers**: the average $\frac1N\sum_{i=1}^NX_i$ of $N$ independent repetitions converges to $E[X]$ as $N\to\infty$ ([Ch. 20](20-statistical-analysis.md)). Using the generalized density of Section 19.7, the two formulas are one: $E[X]=\int xf_X(x)\,dx$ with impulses turning the integral into a sum.

> [!warning] Common Pitfall — the mean need not exist
> The **Cauchy** density $f_X(x)=\dfrac{1}{\pi(1+x^2)}$ — the distribution of $\tan\Theta$ for uniform $\Theta$, and of the ratio of two independent zero-mean Gaussians — is symmetric about 0, yet $\int\lvert x\rvert f_X\,dx=\infty$, so $E[X]$ does not exist. Averages of Cauchy samples do not settle down as $N$ grows: the average of $N$ independent Cauchy variables is again standard Cauchy. Heavy-tailed quantities (ratios of noisy measurements, impulsive interference) can defeat averaging; always check that the mean exists before relying on it.

### Expectation of a function of a random variable

To find $E[g(X)]$ we could first find the distribution of $Y=g(X)$ (Sections 19.7–19.8) and then its mean. Fortunately that is unnecessary.

> [!theorem] Expectation of a function (law of the unconscious statistician)
> $$
> E[g(X)]=\sum_kg(x_k)\,p_X(x_k)\qquad\text{or}\qquad E[g(X)]=\int_{-\infty}^{\infty}g(x)\,f_X(x)\,dx .
> $$
> For two random variables, $E[g(X,Y)]=\iint g(x,y)f_{XY}(x,y)\,dx\,dy$ (or the double sum).

*Proof (discrete case).* $E[Y]=\sum_yy\,P(Y=y)=\sum_yy\sum_{x:g(x)=y}p_X(x)=\sum_xg(x)p_X(x)$, by regrouping the terms of an absolutely convergent sum. The continuous case follows by approximating with discrete variables, or directly from the transformation formula for monotone $g$. $\square$

In Example 19.20, $E[Y]=E[X^2]=4(0.1)+1(0.2)+0+1(0.2)+4(0.1)=1.2$, the same as from the PMF of $Y$: $0(0.4)+1(0.4)+4(0.2)=1.2$. ✓

### Linearity

> [!theorem] Linearity of expectation
> For any random variables $X,Y$ (on the same space, with finite means) and constants $a,b,c$:
> $$
> E[aX+bY+c]=a\,E[X]+b\,E[Y]+c .
> $$
> By induction, $E\left[\sum_iX_i\right]=\sum_iE[X_i]$. **No independence is required.**

*Proof.* By the two-variable form of the previous theorem, $E[aX+bY]=\iint(ax+by)f_{XY}\,dx\,dy=a\int xf_X\,dx+b\int yf_Y\,dy$, using the marginal formulas of Section 19.6. $\square$

Linearity combined with indicator variables is one of the most effective tricks in applied probability. Since $E[I_A]=1\cdot P(A)+0\cdot P(A^c)=P(A)$, the expected number of events that occur among $A_1,\dots,A_n$ is

$$
E\left[\sum_iI_{A_i}\right]=\sum_iP(A_i),
$$

however the events depend on each other.

> [!example] Worked Example 19.38 — Means by indicators
> (a) Find the mean of $X\sim\operatorname{Bin}(n,p)$. (b) In Example 19.5, find the expected number of *pairs* of nodes that share an address, for $n=300$ nodes and $M=65\,536$ addresses.
>
> **(a)** Write $X=I_1+\cdots+I_n$, $I_i$ the indicator of success on trial $i$. Then $E[X]=\sum_iP(\text{success on }i)=np$. (The direct route, $\sum_kk\binom nkp^kq^{n-k}$, needs an index shift and the binomial theorem — and gives the same answer.)
>
> **(b)** For each of the $\binom{300}{2}=44\,850$ pairs, let $I_{ij}$ indicate that nodes $i$ and $j$ collide; $P(I_{ij}=1)=1/M$. The indicators are dependent, but linearity does not care:
> $$
> E[\text{number of colliding pairs}]=\binom{300}{2}\frac{1}{65\,536}=0.684 .
> $$
> Compare: $P(\text{at least one collision})=0.496$. Since $P(N\ge1)\le E[N]$ for any count $N$ (Markov's inequality below), the expected number of collisions is a quick upper bound on the collision probability, and an accurate one when collisions are rare.

### Mean of a non-negative variable and mean time to failure

For a non-negative random variable there is a useful alternative formula in terms of the tail of the distribution.

> [!theorem] Tail formula
> If $T\ge0$, then
> $$
> E[T]=\int_0^\infty P(T>t)\,dt=\int_0^\infty\bigl(1-F_T(t)\bigr)\,dt .
> $$
> For a non-negative integer-valued $N$, $E[N]=\sum_{k=0}^{\infty}P(N>k)$.

*Proof.* Write $T=\int_0^\infty I_{\{T>t\}}\,dt$ (the integrand is 1 for $t<T$ and 0 afterwards) and take expectations inside the integral: $E[T]=\int_0^\infty P(T>t)\,dt$. $\square$

In reliability this is the **mean time to failure** $\text{MTTF}=\int_0^\infty R(t)\,dt$ — the area under the reliability curve.

> [!example] Worked Example 19.39 — MTTF of redundant configurations
> Components have exponential lifetimes with $\lambda=2\times10^{-5}\ \text{h}^{-1}$. Find the MTTF of a single unit, two units in active parallel, TMR (2-out-of-3), and a unit with a perfectly switched cold spare.
>
> **Step 1 — single unit.** $\int_0^\infty e^{-\lambda t}dt=1/\lambda=50\,000$ h.
>
> **Step 2 — active parallel.** $\int_0^\infty(2e^{-\lambda t}-e^{-2\lambda t})\,dt=\frac2\lambda-\frac1{2\lambda}=\frac{3}{2\lambda}=75\,000$ h.
>
> **Step 3 — TMR.** $\int_0^\infty(3e^{-2\lambda t}-2e^{-3\lambda t})\,dt=\frac{3}{2\lambda}-\frac{2}{3\lambda}=\frac{5}{6\lambda}=41\,667$ h.
>
> **Step 4 — cold standby.** $\int_0^\infty(1+\lambda t)e^{-\lambda t}\,dt=\frac1\lambda+\frac1\lambda=\frac2\lambda=100\,000$ h. (Directly: $E[T_1+T_2]=E[T_1]+E[T_2]$ by linearity.)
>
> TMR has a *shorter* mean life than a single unit, although it is more reliable for $\lambda t<\ln2$ (graph in Section 19.4). The MTTF is a poor figure of merit for systems intended for short missions; quote $R(t_{\text{mission}})$ instead.

### Variance and standard deviation

> [!definition] Variance, standard deviation
> The **variance** of $X$ is the mean-square deviation from the mean,
> $$
> \operatorname{Var}(X)=\sigma_X^2=E\bigl[(X-\mu_X)^2\bigr],
> $$
> and the **standard deviation** is $\sigma_X=\sqrt{\operatorname{Var}(X)}$, which has the same units as $X$.

Expanding the square and using linearity gives the computational form

$$
\operatorname{Var}(X)=E[X^2]-2\mu_XE[X]+\mu_X^2=E[X^2]-\mu_X^2 .
$$

> [!theorem] Properties of variance
> 1. $\operatorname{Var}(X)\ge0$, with equality if and only if $X$ is constant with probability 1.
> 2. $\operatorname{Var}(aX+b)=a^2\operatorname{Var}(X)$: a shift does not change spread; a scale multiplies the standard deviation by $\lvert a\rvert$.
> 3. If $X,Y$ are independent, $\operatorname{Var}(X+Y)=\operatorname{Var}(X)+\operatorname{Var}(Y)$; more generally (below) a covariance term appears.

*Proof of 2.* $aX+b-E[aX+b]=a(X-\mu_X)$, so $\operatorname{Var}(aX+b)=E[a^2(X-\mu_X)^2]=a^2\operatorname{Var}(X)$. $\square$

The decomposition $E[X^2]=\mu_X^2+\sigma_X^2$ has a direct electrical meaning. If $X$ is a voltage across $1\ \Omega$, $E[X^2]$ is the average (total) power, $\mu_X^2$ is the power of the DC component and $\sigma_X^2$ is the power of the fluctuating (AC) component; $\sigma_X$ is the AC RMS value read by a true-RMS meter with AC coupling.

> [!example] Worked Example 19.40 — Quantization noise
> An $N$-bit ADC with full-scale range $\pm V_{\text{FS}}$ has step $\Delta=2V_{\text{FS}}/2^N$. For a busy input, the rounding error is modelled as $e\sim\mathcal{U}(-\Delta/2,\Delta/2)$. Find the mean and variance of $e$, and the signal-to-quantization-noise ratio (SQNR) for a full-scale sine wave, evaluated for $N=8$, $V_{\text{FS}}=1$ V.
>
> **Step 1 — mean.** By symmetry $E[e]=0$.
>
> **Step 2 — variance.**
> $$
> \operatorname{Var}(e)=E[e^2]=\frac{1}{\Delta}\int_{-\Delta/2}^{\Delta/2}x^2\,dx=\frac{1}{\Delta}\cdot\frac{2}{3}\left(\frac{\Delta}{2}\right)^3=\frac{\Delta^2}{12}.
> $$
> The RMS quantization noise is $\Delta/\sqrt{12}\approx0.289\Delta$.
>
> **Step 3 — SQNR.** A full-scale sine has power $V_{\text{FS}}^2/2$. With $\Delta=2V_{\text{FS}}/2^N$,
> $$
> \text{SQNR}=\frac{V_{\text{FS}}^2/2}{\Delta^2/12}=\frac{3}{2}\,2^{2N},\qquad 10\log_{10}\text{SQNR}=6.02N+1.76\ \text{dB}.
> $$
>
> **Step 4 — numbers.** $\Delta=2/256=7.8125$ mV, $\operatorname{Var}(e)=5.086\times10^{-6}\ \text{V}^2$, RMS noise $2.255$ mV, SQNR $=10\log_{10}(0.5/5.086\times10^{-6})=49.93$ dB, matching $6.02(8)+1.76=49.92$ dB. Every extra bit buys about 6 dB.

### Means and variances of the key distributions

| Distribution | Mean | Variance | Derivation hint |
|---|---|---|---|
| Bernoulli($p$) | $p$ | $p(1-p)$ | $E[X^2]=E[X]=p$ |
| $\operatorname{Bin}(n,p)$ | $np$ | $np(1-p)$ | sum of $n$ independent Bernoullis |
| $\operatorname{Geom}(p)$ | $1/p$ | $(1-p)/p^2$ | tail formula: $\sum_k(1-p)^k=1/p$ |
| $\operatorname{Poisson}(\lambda)$ | $\lambda$ | $\lambda$ | MGF below, or $E[X(X-1)]=\lambda^2$ |
| $\mathcal{U}(a,b)$ | $(a+b)/2$ | $(b-a)^2/12$ | Example 19.40 shifted and scaled |
| $\operatorname{Exp}(\lambda)$ | $1/\lambda$ | $1/\lambda^2$ | $E[T^n]=n!/\lambda^n$ |
| $\mathcal{N}(\mu,\sigma^2)$ | $\mu$ | $\sigma^2$ | standardize; $E[Z^2]=1$ by parts |
| Rayleigh($\sigma$) | $\sigma\sqrt{\pi/2}$ | $(2-\pi/2)\sigma^2$ | $E[R^2]=2\sigma^2$ |

*Two derivations.* For $Z\sim\mathcal{N}(0,1)$, $E[Z]=0$ by symmetry and, integrating by parts with $u=z$, $dv=z\varphi(z)\,dz$ (so $v=-\varphi(z)$): $E[Z^2]=\int z^2\varphi(z)\,dz=[-z\varphi(z)]_{-\infty}^{\infty}+\int\varphi(z)\,dz=1$. Then $X=\mu+\sigma Z$ has mean $\mu$ and variance $\sigma^2$, which justifies the names of the parameters. For $T\sim\operatorname{Exp}(\lambda)$, repeated integration by parts gives $E[T^n]=\int_0^\infty t^n\lambda e^{-\lambda t}dt=n!/\lambda^n$, so $E[T]=1/\lambda$, $E[T^2]=2/\lambda^2$ and $\operatorname{Var}(T)=1/\lambda^2$: the standard deviation of an exponential lifetime equals its mean.

### Sums, covariance and independence

For two random variables, expanding $E[(X+Y-\mu_X-\mu_Y)^2]$ gives

$$
\operatorname{Var}(X+Y)=\operatorname{Var}(X)+\operatorname{Var}(Y)+2\operatorname{Cov}(X,Y),\qquad \operatorname{Cov}(X,Y)=E[(X-\mu_X)(Y-\mu_Y)]=E[XY]-\mu_X\mu_Y .
$$

> [!theorem] Expectation of a product of independent RVs
> If $X$ and $Y$ are independent, then $E[XY]=E[X]E[Y]$, more generally $E[g(X)h(Y)]=E[g(X)]E[h(Y)]$. Hence $\operatorname{Cov}(X,Y)=0$ and variances of independent RVs add:
> $$
> \operatorname{Var}\left(\sum_{i=1}^nX_i\right)=\sum_{i=1}^n\operatorname{Var}(X_i)\qquad(X_i\text{ independent}).
> $$

*Proof.* $E[g(X)h(Y)]=\iint g(x)h(y)f_X(x)f_Y(y)\,dx\,dy=\int gf_X\,dx\int hf_Y\,dy$, the double integral factoring because the joint density does. $\square$

The converse is false: zero covariance (**uncorrelated**) does not imply independence. For $X\sim\mathcal{N}(0,1)$ and $Y=X^2$, $\operatorname{Cov}(X,Y)=E[X^3]=0$, yet $Y$ is a function of $X$. Covariance, the correlation coefficient and their estimation from data are developed in [Ch. 20](20-statistical-analysis.md); autocorrelation of signals in [Ch. 21](21-random-processes.md).

Two consequences are used constantly. First, the **variance of an average**: if $X_1,\dots,X_N$ are independent with common variance $\sigma^2$, then $\bar X=\frac1N\sum X_i$ has $\operatorname{Var}(\bar X)=\sigma^2/N$, so averaging $N$ independent readings reduces the RMS noise by $\sqrt N$ (signal averaging gives 3 dB of SNR per doubling of $N$). Second, **tolerance stack-up**: independent tolerances add in quadrature, not linearly.

> [!example] Worked Example 19.41 — DC plus noise, and resistors in series
> (a) A sensor output is $V=2+N$ V with $N\sim\mathcal{N}(0,0.5^2)$. Find the mean, variance and mean-square value of $V$, and the average power it delivers to $1\ \Omega$. (b) Two independent resistors, each with mean $1000\ \Omega$ and standard deviation $10\ \Omega$, are connected in series. Find the mean and standard deviation of the total, absolute and relative.
>
> **(a)** $E[V]=2$ V, $\operatorname{Var}(V)=0.25\ \text{V}^2$, $E[V^2]=2^2+0.25=4.25\ \text{V}^2$. Power into $1\ \Omega$: $4.25$ W, of which $4$ W is DC and $0.25$ W is noise.
>
> **(b)** Mean $2000\ \Omega$; variance $10^2+10^2=200\ \Omega^2$; standard deviation $\sqrt{200}=14.14\ \Omega$. Relative spread: $14.14/2000=0.71\,\%$, compared with $1\,\%$ for each resistor. Independent errors partially cancel; a worst-case (linear) stack-up of $\pm20\ \Omega$ would be pessimistic.

### Conditional expectation and the law of total expectation

Conditioning works for expectations exactly as for probabilities. For an event $B$, $E[X\mid B]$ is the mean computed from the conditional distribution of $X$ given $B$. If $B_1,\dots,B_n$ partition $\Omega$,

$$
E[X]=\sum_{i=1}^nE[X\mid B_i]\,P(B_i),
$$

the **law of total expectation**, proved by applying total probability inside the sum or integral that defines $E[X]$.

> [!example] Worked Example 19.42 — Mean number of transmissions by conditioning
> In the ARQ protocol of Example 19.23, each attempt succeeds with probability $p$. Find $E[N]$ without summing a series.
>
> **Step 1 — condition on the first attempt.** If it succeeds (probability $p$), $N=1$. If it fails (probability $1-p$), one attempt has been used and, by independence, the situation restarts: $E[N\mid\text{fail}]=1+E[N]$.
>
> **Step 2 — solve.** $E[N]=p\cdot1+(1-p)(1+E[N])$, so $pE[N]=1$ and $E[N]=1/p$.
>
> **Step 3 — numbers.** With $p=0.8$, $E[N]=1.25$ attempts: the throughput of the link is $80\,\%$ of its raw rate. Check by the tail formula: $\sum_{k\ge0}(0.2)^k=1/0.8=1.25$. ✓

### Inequalities: Markov and Chebyshev

Means and variances bound tail probabilities even when the distribution is unknown.

> [!theorem] Markov and Chebyshev inequalities
> **Markov:** if $X\ge0$ and $a>0$, then $P(X\ge a)\le\dfrac{E[X]}{a}$.
>
> **Chebyshev:** for any $X$ with mean $\mu$ and variance $\sigma^2$, and any $k>0$,
> $$
> P(\lvert X-\mu\rvert\ge k\sigma)\le\frac{1}{k^2}.
> $$

*Proof.* Markov: $X\ge a\,I_{\{X\ge a\}}$ for $X\ge0$; take expectations, $E[X]\ge aP(X\ge a)$. Chebyshev: apply Markov to the non-negative variable $(X-\mu)^2$ with $a=k^2\sigma^2$. $\square$

> [!example] Worked Example 19.43 — Guaranteed yield without a distribution
> A supplier guarantees only that a batch of resistors has mean $1000\ \Omega$ and standard deviation $10\ \Omega$. What fraction can lie outside $1000\pm30\ \Omega$? Compare with the Gaussian model.
>
> **Step 1 — Chebyshev with $k=3$.** $P(\lvert R-1000\rvert\ge30)\le1/9=0.111$.
>
> **Step 2 — Gaussian.** $2Q(3)=0.0027$.
>
> Chebyshev's bound is distribution-free and therefore very conservative — about 40 times the Gaussian value here — but it is a *guarantee*. It is attained only by a three-point distribution, which real batches rarely resemble. Chebyshev's inequality also proves the weak law of large numbers in one line: $P(\lvert\bar X-\mu\rvert\ge\epsilon)\le\sigma^2/(N\epsilon^2)\to0$.

```desmos-graph
left=1; right=5; bottom=0; top=1.05
---
y=\frac{1}{x^{2}}|red|1<x<5
y=1-\operatorname{erf}\left(\frac{x}{\sqrt{2}}\right)|blue|1<x<5
```

*Chebyshev bound $1/k^2$ (red) on $P(\lvert X-\mu\rvert\ge k\sigma)$, and the exact value $2Q(k)$ for a Gaussian (blue). The distribution-free bound decays only polynomially; Gaussian tails decay like $e^{-k^2/2}$.*

### Moments and the shape of a distribution

> [!definition] Moments
> The **$n$-th moment** of $X$ is $m_n=E[X^n]$ and the **$n$-th central moment** is $\mu_n=E[(X-\mu)^n]$. Thus $m_1=\mu$, $\mu_2=\sigma^2$. The normalized third and fourth central moments are the **skewness** and **kurtosis**:
> $$
> \gamma_1=\frac{\mu_3}{\sigma^3},\qquad \beta_2=\frac{\mu_4}{\sigma^4}\quad(\text{excess kurtosis } \beta_2-3).
> $$

Skewness measures asymmetry (positive when the right tail is longer); kurtosis measures the weight of the tails relative to a Gaussian, for which $\gamma_1=0$ and $\beta_2=3$. Central moments are found from ordinary moments by expanding the binomial: for instance $\mu_3=m_3-3\mu m_2+2\mu^3$. For the exponential distribution, with $m_n=n!/\lambda^n$, this gives $\mu_3=(6-6+2)/\lambda^3=2/\lambda^3$ and $\gamma_1=2$; similarly $\beta_2=9$ (excess kurtosis 6). Impulsive noise in power-line communications and switching transients is characterized by large excess kurtosis — its tails are much heavier than Gaussian, so Gaussian-based error-rate predictions are optimistic.

Knowing all the moments usually determines the distribution (it does, for example, whenever the moment-generating function below exists near $0$). The transforms in the rest of the section encode the whole sequence of moments in one function.

### The moment-generating function

> [!definition] Moment-generating function (MGF)
> $$
> M_X(s)=E\bigl[e^{sX}\bigr]=\int_{-\infty}^{\infty}e^{sx}f_X(x)\,dx\quad\Bigl(\text{or }\sum_ke^{sx_k}p_X(x_k)\Bigr),
> $$
> defined for the real $s$ for which the expectation is finite. It always exists at $s=0$, where $M_X(0)=1$.

The MGF of a continuous RV is the two-sided Laplace transform ([Ch. 14](../part-5-transforms/14-laplace-transforms.md)) of $f_X$ evaluated at $-s$. Its name comes from the following property.

> [!theorem] Properties of the MGF
> If $M_X(s)$ is finite in an open interval around $s=0$, then:
> 1. **Moments:** $E[X^n]=M_X^{(n)}(0)=\dfrac{d^nM_X}{ds^n}\Big|_{s=0}$, and $M_X(s)=\sum_{n=0}^{\infty}\dfrac{E[X^n]}{n!}s^n$.
> 2. **Uniqueness:** $M_X$ determines the distribution of $X$.
> 3. **Linear transformation:** $M_{aX+b}(s)=e^{bs}M_X(as)$.
> 4. **Sums:** if $X,Y$ are independent, $M_{X+Y}(s)=M_X(s)M_Y(s)$.

*Proof of 1 and 4.* Expand $e^{sX}=\sum_n(sX)^n/n!$ and take expectations term by term (justified by the finiteness of $M_X$ near 0); the coefficient of $s^n/n!$ is $E[X^n]$. For 4, $E[e^{s(X+Y)}]=E[e^{sX}e^{sY}]=E[e^{sX}]E[e^{sY}]$ by the product theorem for independent variables. $\square$

A related convenience is the **cumulant** function $\ln M_X(s)$: its first two derivatives at 0 are the mean and the variance directly, $(\ln M_X)'(0)=\mu$, $(\ln M_X)''(0)=\sigma^2$.

> [!example] Worked Example 19.44 — MGFs of the exponential, Poisson and Gaussian
> Find the MGFs and use them to obtain the means and variances.
>
> **Exponential.** For $s<\lambda$: $M_T(s)=\int_0^\infty e^{st}\lambda e^{-\lambda t}dt=\dfrac{\lambda}{\lambda-s}$. Then $M'(s)=\lambda/(\lambda-s)^2$, $M''(s)=2\lambda/(\lambda-s)^3$, so $E[T]=1/\lambda$, $E[T^2]=2/\lambda^2$, $\operatorname{Var}(T)=1/\lambda^2$. ✓ (The MGF does not exist for $s\ge\lambda$ — the exponential tail is not light enough.)
>
> **Poisson.** $M_X(s)=\sum_ke^{sk}e^{-\lambda}\lambda^k/k!=e^{-\lambda}\sum_k(\lambda e^s)^k/k!=\exp\bigl(\lambda(e^s-1)\bigr)$. Its logarithm is $\lambda(e^s-1)$, whose first and second derivatives at 0 are both $\lambda$: **mean = variance = $\lambda$**.
>
> **Gaussian.** For $Z\sim\mathcal{N}(0,1)$, complete the square: $sz-z^2/2=-\frac12(z-s)^2+\frac{s^2}{2}$, so $M_Z(s)=e^{s^2/2}\int\varphi(z-s)\,dz=e^{s^2/2}$. By property 3, $X=\mu+\sigma Z$ has
> $$
> M_X(s)=\exp\!\left(\mu s+\tfrac12\sigma^2s^2\right).
> $$
> $\ln M_X=\mu s+\frac12\sigma^2s^2$, giving mean $\mu$ and variance $\sigma^2$ at once. Expanding $e^{s^2/2}=\sum s^{2k}/(2^kk!)$ also gives all the moments of $Z$: odd moments vanish and $E[Z^{2k}]=(2k)!/(2^kk!)=1\cdot3\cdot5\cdots(2k-1)$, so $E[Z^4]=3$ (kurtosis 3).

> [!example] Worked Example 19.45 — Closure of Poisson and Gaussian families under addition
> (a) Packets from two independent sources arrive as $\operatorname{Poisson}(\lambda_1)$ and $\operatorname{Poisson}(\lambda_2)$ counts per slot. What is the distribution of the total? (b) A receiver input has independent thermal noise $N_1\sim\mathcal{N}(0,1)$ and interference $N_2\sim\mathcal{N}(0,4)$ (in $\mu\text{V}^2$). What is the distribution of $N_1+N_2$, and the probability that it exceeds $5\ \mu$V?
>
> **(a)** $M(s)=e^{\lambda_1(e^s-1)}e^{\lambda_2(e^s-1)}=e^{(\lambda_1+\lambda_2)(e^s-1)}$, the MGF of $\operatorname{Poisson}(\lambda_1+\lambda_2)$. By uniqueness, merged independent Poisson streams are Poisson with the summed rate. (Doing this by discrete convolution requires the binomial theorem; the MGF makes it a one-liner.)
>
> **(b)** $M(s)=e^{s^2/2}e^{4s^2/2}=e^{5s^2/2}$, so $N_1+N_2\sim\mathcal{N}(0,5)$. Then $P(N_1+N_2>5)=Q(5/\sqrt5)=Q(2.236)=0.0127$. Independent Gaussian noises add in *power* (variance), and the sum stays Gaussian — the property that makes Gaussian noise analysis of linear systems so clean ([Ch. 21](21-random-processes.md)).

### The characteristic function

The MGF may fail to exist (for the Cauchy distribution it is infinite for every $s\ne0$). Replacing the real $s$ by an imaginary $j\omega$ cures this, because $\lvert e^{j\omega X}\rvert=1$.

> [!definition] Characteristic function
> $$
> \Phi_X(\omega)=E\bigl[e^{j\omega X}\bigr]=\int_{-\infty}^{\infty}f_X(x)\,e^{j\omega x}\,dx ,\qquad\omega\in\mathbb{R}.
> $$
> It exists for every random variable, with $\lvert\Phi_X(\omega)\rvert\le\Phi_X(0)=1$.

Comparing with the Fourier transform of [Ch. 16](../part-5-transforms/16-fourier-transforms.md), $\Phi_X(\omega)$ is the Fourier transform of the density with the sign of the exponent reversed: $\Phi_X(\omega)=\mathcal{F}\{f_X\}(-\omega)$. Every Fourier result transfers immediately:

- **Inversion:** $f_X(x)=\dfrac{1}{2\pi}\displaystyle\int_{-\infty}^{\infty}\Phi_X(\omega)e^{-j\omega x}\,d\omega$, so $\Phi_X$ determines the distribution.
- **Convolution theorem:** for independent $X,Y$, $\Phi_{X+Y}(\omega)=\Phi_X(\omega)\Phi_Y(\omega)$ — the transform-domain version of $f_{X+Y}=f_X*f_Y$ (Section 19.6).
- **Moments:** if $E[\lvert X\rvert^n]<\infty$, $E[X^n]=j^{-n}\Phi_X^{(n)}(0)$; smoothness of $\Phi_X$ at the origin reflects how many moments exist, just as smoothness of a signal reflects the decay of its spectrum.
- **Symmetry:** $\Phi_X(-\omega)=\Phi_X(\omega)^*$; $\Phi_X$ is real and even exactly when $f_X$ is even.

Standard pairs: $\mathcal{N}(\mu,\sigma^2)\mapsto e^{j\mu\omega-\sigma^2\omega^2/2}$; $\operatorname{Exp}(\lambda)\mapsto\lambda/(\lambda-j\omega)$; $\mathcal{U}(-a,a)\mapsto\sin(a\omega)/(a\omega)$; Cauchy $\mapsto e^{-\lvert\omega\rvert}$ (not differentiable at 0, matching the non-existent mean); $\operatorname{Poisson}(\lambda)\mapsto\exp(\lambda(e^{j\omega}-1))$.

For integer-valued random variables one also uses the **probability generating function** $G_N(z)=E[z^N]=\sum_kp_N(k)z^k$, which is the Z-transform of the PMF ([Ch. 17](../part-5-transforms/17-z-transforms.md)) evaluated at $z^{-1}$; $G_N'(1)=E[N]$ and the PGF of a sum of independent counts is the product of the PGFs.

> [!example] Worked Example 19.46 — Sum of independent uniform errors via characteristic functions
> A measurement accumulates $n$ independent rounding errors $e_i\sim\mathcal{U}(-\frac12,\frac12)$ (in units of one LSB). Find the characteristic function, mean and variance of the total error $S_n=e_1+\cdots+e_n$, and compare its distribution with a Gaussian for $n=3$.
>
> **Step 1 — one error.** $\Phi_e(\omega)=\int_{-1/2}^{1/2}e^{j\omega x}dx=\dfrac{\sin(\omega/2)}{\omega/2}$.
>
> **Step 2 — the sum.** By independence, $\Phi_{S_n}(\omega)=\left(\dfrac{\sin(\omega/2)}{\omega/2}\right)^n$ — the $n$-fold product, corresponding to the $n$-fold convolution of rectangles (the triangle of Example 19.27 for $n=2$).
>
> **Step 3 — moments.** Expand near $\omega=0$: $\frac{\sin u}{u}=1-\frac{u^2}{6}+\cdots$ with $u=\omega/2$, so $\Phi_e(\omega)=1-\frac{\omega^2}{24}+\cdots$. Matching with $\Phi(\omega)=1+j\omega E[e]-\frac{\omega^2}{2}E[e^2]+\cdots$ gives $E[e]=0$, $E[e^2]=1/12$ (as in Example 19.40). Then $\Phi_{S_n}(\omega)=1-\frac{n\omega^2}{24}+\cdots$, so $\operatorname{Var}(S_n)=n/12$.
>
> **Step 4 — Gaussian comparison.** For $n=3$, $\operatorname{Var}(S_3)=1/4$, $\sigma=1/2$. The exact density (three convolved rectangles, a piecewise quadratic) has $f_{S_3}(0)=3/4$; the Gaussian $\mathcal{N}(0,1/4)$ has $f(0)=1/(\sqrt{2\pi}\cdot\frac12)=0.798$ — already within $6\,\%$. Since $\left(\frac{\sin(\omega/2)}{\omega/2}\right)^n\approx\left(1-\frac{\omega^2}{24}\right)^n\to$ a Gaussian shape after rescaling by $\sqrt n$, the characteristic function makes the central limit theorem plausible in a single line; [Ch. 20](20-statistical-analysis.md) states it precisely.

```desmos-graph
left=-2.2; right=2.2; bottom=0; top=0.9
---
y=\frac{1}{2}(x+1.5)^{2}|blue|-1.5<x<-0.5
y=0.75-x^{2}|blue|-0.5<x<0.5
y=\frac{1}{2}(1.5-x)^{2}|blue|0.5<x<1.5
y=\frac{1}{0.5\sqrt{2\pi}}e^{-\frac{x^{2}}{0.5}}|red|dashed
```

*Density of the sum of three independent $\mathcal{U}(-\frac12,\frac12)$ errors (blue, piecewise quadratic) and the Gaussian with the same variance $1/4$ (red, dashed). Three uniform terms are already nearly Gaussian; the difference lies mainly in the tails, which the uniform sum truncates at $\pm1.5$.*

> [!warning] Common Pitfall — sign conventions for transforms
> The characteristic function uses $e^{+j\omega x}$, the Fourier transform of this book uses $e^{-j\omega t}$, and the MGF uses the real exponent $e^{sx}$. Mixing conventions flips the sign of odd moments (and of the phase). When borrowing a Fourier pair for a characteristic function, replace $\omega$ by $-\omega$; for even densities nothing changes.

> [!tip] Engineering Note — moments as the designer's summary
> Mean (DC level, bias), variance (AC power, noise power, tolerance), and the shape parameters (skewness, kurtosis) are the quantities that appear in specifications and link budgets. Second-order reasoning — means and variances, which add for independent contributions — handles the bulk of noise and tolerance analysis; the full distribution, usually Gaussian by the central limit theorem, is needed only at the final step, when the variance is converted into an error probability with $Q(\cdot)$.

## Summary

| Concept | Key formula |
|---|---|
| Axioms | $P(A)\ge0$, $P(\Omega)=1$, $P(\bigcup A_i)=\sum P(A_i)$ for disjoint $A_i$ |
| Complement, addition rule | $P(A^c)=1-P(A)$, $P(A\cup B)=P(A)+P(B)-P(AB)$ |
| Union bound | $P(\bigcup A_i)\le\sum P(A_i)$ |
| Counting | $n^k$, $n!/(n-k)!$, $\binom nk=\dfrac{n!}{k!(n-k)!}$, $\dbinom{n+k-1}{k}$ |
| Conditional probability | $P(A\mid B)=P(AB)/P(B)$ |
| Chain rule | $P(A_1\cdots A_n)=P(A_1)P(A_2\mid A_1)\cdots P(A_n\mid A_1\cdots A_{n-1})$ |
| Total probability | $P(A)=\sum_iP(A\mid B_i)P(B_i)$ |
| Bayes' theorem | $P(B_k\mid A)=P(A\mid B_k)P(B_k)/\sum_iP(A\mid B_i)P(B_i)$ |
| Odds form | posterior odds $=$ likelihood ratio $\times$ prior odds |
| MAP threshold, $\pm A$ in $\mathcal{N}(0,\sigma^2)$ | $x_{\text{th}}=\dfrac{\sigma^2}{2A}\ln\dfrac{p_0}{p_1}$ |
| Independence | $P(AB)=P(A)P(B)$; for RVs $f_{XY}=f_Xf_Y$ |
| Series / parallel reliability | $\prod p_i$ / $1-\prod(1-p_i)$ |
| $k$-out-of-$n$, TMR | $\sum_{i\ge k}\binom ni p^i(1-p)^{n-i}$; $3p^2-2p^3$ |
| Binomial, geometric, Poisson | $\binom nkp^kq^{n-k}$; $q^{k-1}p$; $e^{-\lambda}\lambda^k/k!$ |
| Uniform, exponential | $1/(b-a)$; $\lambda e^{-\lambda t}$, $P(T>t)=e^{-\lambda t}$ |
| Gaussian | $\dfrac{1}{\sqrt{2\pi}\sigma}e^{-(x-\mu)^2/2\sigma^2}$, $P(X>x)=Q\bigl((x-\mu)/\sigma\bigr)$ |
| Rayleigh | $\dfrac{r}{\sigma^2}e^{-r^2/2\sigma^2}$, $P(R>r)=e^{-r^2/2\sigma^2}$ |
| CDF | $F_X(x)=P(X\le x)$, $P(a<X\le b)=F_X(b)-F_X(a)$, $f_X=F_X'$ |
| Marginal, conditional | $f_X=\int f_{XY}\,dy$, $f_{Y\mid X}=f_{XY}/f_X$ |
| Sum of independent RVs | $f_{X+Y}=f_X*f_Y$ |
| Density transformation | $f_Y(y)=\sum_if_X(x_i)/\lvert g'(x_i)\rvert$, $g(x_i)=y$ |
| Max / min of independent RVs | $F_{\max}=\prod F_i$, $P(\min>x)=\prod(1-F_i)$ |
| Hazard rate | $h(t)=f_T(t)/R(t)$, $R(t)=\exp\left(-\int_0^th\right)$ |
| Expectation, LOTUS | $E[g(X)]=\int g(x)f_X(x)\,dx$; $E[aX+bY]=aE[X]+bE[Y]$ |
| Tail formula, MTTF | $E[T]=\int_0^\infty P(T>t)\,dt$ |
| Variance | $\operatorname{Var}(X)=E[X^2]-\mu^2$, $\operatorname{Var}(aX+b)=a^2\operatorname{Var}(X)$ |
| Sum of variances | $\operatorname{Var}(X+Y)=\operatorname{Var}X+\operatorname{Var}Y+2\operatorname{Cov}(X,Y)$ |
| Quantization noise | $\Delta^2/12$; SQNR $=6.02N+1.76$ dB |
| Markov, Chebyshev | $P(X\ge a)\le E[X]/a$; $P(\lvert X-\mu\rvert\ge k\sigma)\le1/k^2$ |
| MGF | $M_X(s)=E[e^{sX}]$, $E[X^n]=M_X^{(n)}(0)$, $M_{X+Y}=M_XM_Y$ |
| Characteristic function | $\Phi_X(\omega)=E[e^{j\omega X}]$, $\Phi_{X+Y}=\Phi_X\Phi_Y$ |
| Gaussian MGF | $\exp(\mu s+\sigma^2s^2/2)$ |

## Further reading

- A. Leon-Garcia, *Probability, Statistics, and Random Processes for Electrical Engineering*, 3rd ed., Pearson — the closest match to this chapter, with many communications and reliability examples.
- A. Papoulis and S. U. Pillai, *Probability, Random Variables, and Stochastic Processes*, 4th ed., McGraw-Hill — the standard graduate reference, especially for functions of random variables and characteristic functions.
- D. P. Bertsekas and J. N. Tsitsiklis, *Introduction to Probability*, 2nd ed., Athena Scientific — exceptionally clear on conditioning, Bayes' rule and expectation.
- S. M. Ross, *A First Course in Probability*, 10th ed., Pearson — combinatorics and classic problems in depth.
