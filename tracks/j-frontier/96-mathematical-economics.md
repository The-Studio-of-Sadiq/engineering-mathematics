# Chapter 96 — Mathematical Economics and Social Choice

[← Previous](95-representation-theory.md) · [Track J](README.md) · [Table of Contents](../../00-front/table-of-contents.md) · [Next →](97-measure-theory.md)

> **Level 3 · Track J — Mathematics at the Frontier**
> **Builds on:** [Ch. 6](../../part-2-calculus/06-multivariable-calculus.md), [Ch. 23](../../part-7-numerical-computational/23-optimization.md), [Ch. 42](../c-game-theory/42-mechanism-design.md), [Ch. 44](../c-game-theory/44-coalitions-and-voting.md) · **Leads to:** [Ch. 97](97-measure-theory.md), [Ch. 98](98-mathematical-research.md) · **Related:** [Ch. 39](../c-game-theory/39-game-theory-basics.md), [Ch. 41](../c-game-theory/41-incomplete-information-games.md), [Ch. 43](../c-game-theory/43-auctions-and-markets.md), [Ch. 47](../d-optimization/47-convexity.md), [Ch. 30](../a-algorithms/30-limits-of-computation.md)
> **Reading time:** ~30 min

## 96.1 The hook

In 1951 a young economist named Kenneth Arrow published two results that seem to contradict each other.

The first was a mathematical version of Adam Smith's "invisible hand". Let millions of selfish people trade at market prices, each caring only about their own basket of goods. If prices settle where supply meets demand, the outcome is **Pareto efficient**: no central planner, however wise, could rearrange the goods to make one person better off without making someone else worse off. Selfish choices, aggregated by prices, give a collectively optimal result.

The second was his doctoral thesis, *Social Choice and Individual Values*. It asked whether a society could aggregate its members' preferences by *voting*, and it said no. Any rule that turns individual rankings of three or more options into a social ranking, respects unanimity, and judges each pair of options by individual opinions on that pair alone, must be a **dictatorship**: one voter's ranking decides everything.

So markets aggregate preferences perfectly and ballots cannot aggregate them at all. Arrow shared the 1972 Nobel prize largely for these two results, and it is not obvious that they can both be true. Hasn't a market just counted votes in dollars?

This chapter proves both theorems, and on the way it settles several debts left by the game-theory track: why Gibbard–Satterthwaite forces a dictator ([Ch. 42](../c-game-theory/42-mechanism-design.md)), why no trading rule can be efficient and balanced at once, and why haggling ends at the Nash bargaining solution ([Ch. 44](../c-game-theory/44-coalitions-and-voting.md)). The answer to the hook turns out to be about *information*. A price carries more than a vote does, and it asks for less.

## 96.2 The big idea

Start with one person. Economics doesn't assume we can measure happiness. It assumes only that you can *compare*: given two bundles, you prefer one, or the other, or you're indifferent. A utility function is merely a way of writing those comparisons as numbers, and any increasing relabelling of it ($u\mapsto u^3$, $u\mapsto\ln u$) describes the same person. Utility is **ordinal**.

Now two people, Ann and Bob, and two goods, apples and bread. Draw a rectangle whose width is the total number of apples and whose height is the total bread. A point inside is a way to split everything: measure Ann's holdings from the bottom-left corner and Bob's from the top-right. This is the **Edgeworth box**, and every question about two-person exchange becomes geometry inside it.

Through the starting point (the **endowment**) runs one of Ann's indifference curves and one of Bob's. Between them is a lens of allocations that *both* prefer. Trading inside the lens helps both. Trading stops being mutually helpful only when the two indifference curves touch without crossing, that is, when their slopes (the rates at which each would swap apples for bread) are equal. The locus of all such tangencies is the **contract curve**, the set of Pareto-efficient splits.

A price line through the endowment is a budget: at price ratio $p$, Ann can sell apples for bread along it. Each person picks their favourite point on the line, so each sets their own rate of substitution equal to the *same* price. That is why markets hit the contract curve. Everyone's private trade-off is aligned to one public number, and nobody needs to know anyone else's preferences.

Voting fails to do this because a ranking carries no rate of exchange. Arrow's independence axiom says society's verdict on $a$ versus $b$ may depend only on who prefers $a$ to $b$, not on *how much*. Prices express intensity, and ballots by construction don't. Both theorems are about how much information an aggregation rule is allowed to use.

## 96.3 The mathematics

### Preferences and utility

> [!definition] Preference relation
> A **preference** on a set $X$ of bundles is a relation $\succsim$ ("at least as good as") that is **complete** (for all $x,y$: $x\succsim y$ or $y\succsim x$) and **transitive**. A function $u:X\to\mathbb{R}$ **represents** it if $x\succsim y\iff u(x)\ge u(y)$.

On a finite set every preference has a utility (rank the bundles and number them). On $\mathbb{R}^2$ the natural-looking **lexicographic** order, "more apples is better, and bread only breaks ties", has none.

*Proof.* Suppose $u$ represents it. For each apple quantity $x\in\mathbb{R}$, $(x,1)\succ(x,0)$, so the interval $I_x=(u(x,0),u(x,1))$ is non-empty and contains a rational $q_x$. If $x<x'$ then $(x,1)\prec(x',0)$, so $I_x$ lies entirely below $I_{x'}$ and $q_x\neq q_{x'}$. That gives an injection from $\mathbb{R}$ into $\mathbb{Q}$, which is impossible because $\mathbb{R}$ is uncountable ([Ch. 64](../f-logic-foundations/64-set-theory-and-infinity.md)). ∎

Debreu proved in 1954 that **continuity** (strict-preference sets are open) is enough to guarantee a continuous utility, which is what licenses calculus in economics.

### Exchange, efficiency and equilibrium

Let consumers $i=1,\dots,I$ have endowments $\boldsymbol\omega_i\in\mathbb{R}^L_{\ge0}$ of $L$ goods. An **allocation** $(\mathbf{x}_i)$ is **feasible** if $\sum_i\mathbf{x}_i=\sum_i\boldsymbol\omega_i$. It is **Pareto efficient** if no feasible allocation makes everyone weakly better off and someone strictly better off, the multi-objective optimality of [§23.9](../../part-7-numerical-computational/23-optimization.md#239-multi-objective-optimization).

> [!definition] Walrasian (competitive) equilibrium
> A price vector $\mathbf{p}$ and a feasible allocation $(\mathbf{x}^*_i)$ such that each $\mathbf{x}^*_i$ is a best bundle for $i$ in the budget set $\{\mathbf{x}\ge0:\mathbf{p}\cdot\mathbf{x}\le\mathbf{p}\cdot\boldsymbol\omega_i\}$.

> [!example] Ann and Bob
> Ten apples and ten loaves in total. Ann starts with $(8,2)$ and has $u_A=x^{0.6}y^{0.4}$; Bob starts with $(2,8)$ and has $u_B=x^{0.4}y^{0.6}$. Fix the bread price at $1$ and the apple price at $p$.
> Maximising $x^{a}y^{1-a}$ on a budget $I$ by a Lagrange multiplier ([§6.11](../../part-2-calculus/06-multivariable-calculus.md#611-lagrange-multipliers)) gives the Cobb–Douglas rule "spend the fraction $a$ on apples": $x=aI/p$. Ann's income is $8p+2$ and Bob's is $2p+8$. Apples clear when
> $$
> \frac{0.6(8p+2)+0.4(2p+8)}{p}=10\iff 5.6p+4.4=10p\iff p=1 .
> $$
> So Ann ends with $(6,4)$ and Bob with $(4,6)$. Equal rates of substitution, $\frac{0.6}{0.4}\frac{y}{x}=\frac{0.4}{0.6}\frac{10-y}{10-x}$, give the contract curve $y=8x/(18-x)$, and $(6,4)$ lies on it: $48/12=4$.

![An Edgeworth box with Ann's origin at bottom left and Bob's at top right; blue and red indifference curves through the endowment (8,2) enclose a green lens; a purple contract curve runs from corner to corner; the dashed budget line x+y=10 touches both indifference curves at the equilibrium star (6,4)](../../images/ch96-edgeworth-box.svg)

*Figure 96.1 — The Edgeworth box. Trades in the green lens help both. The price line through the endowment meets the contract curve exactly where both indifference curves are tangent to it.*

> [!theorem] First fundamental theorem of welfare economics
> If preferences are **locally non-satiated** (near any bundle there is a strictly better one), every Walrasian equilibrium allocation is Pareto efficient.

*Proof.* Let $(\mathbf{p},\mathbf{x}^*)$ be an equilibrium and suppose a feasible $(\mathbf{y}_i)$ Pareto-dominates it. If $\mathbf{y}_j\succ_j\mathbf{x}^*_j$, then $\mathbf{y}_j$ was unaffordable, so $\mathbf{p}\cdot\mathbf{y}_j>\mathbf{p}\cdot\boldsymbol\omega_j$. If $\mathbf{y}_i\succsim_i\mathbf{x}^*_i$, then $\mathbf{p}\cdot\mathbf{y}_i\ge\mathbf{p}\cdot\boldsymbol\omega_i$. Otherwise $\mathbf{y}_i$ would cost strictly less than the budget, a slightly better bundle nearby would still be affordable, and it would beat $\mathbf{x}^*_i$, contradicting optimality. Summing over all consumers, with at least one strict inequality,

$$
\mathbf{p}\cdot\sum_i\mathbf{y}_i>\mathbf{p}\cdot\sum_i\boldsymbol\omega_i ,
$$

which contradicts feasibility $\sum_i\mathbf{y}_i=\sum_i\boldsymbol\omega_i$. ∎

Notice what the proof did *not* use: convexity, differentiability, or anyone knowing anyone else's tastes. Prices did all the coordinating. Two companions complete the picture. The **second welfare theorem** says that with convex preferences *every* Pareto optimum is an equilibrium after suitable lump-sum transfers. Its proof is a separating hyperplane, whose normal vector *is* the price ([Ch. 47](../d-optimization/47-convexity.md), [Ch. 48](../d-optimization/48-duality.md)). **Existence** of an equilibrium (Arrow and Debreu, 1954) is a fixed-point theorem: the map "raise prices of over-demanded goods" sends the price simplex continuously into itself, the same kind of fixed-point argument (Brouwer's theorem, or Kakutani's extension of it) that gives Nash equilibria ([Ch. 39](../c-game-theory/39-game-theory-basics.md)).

![Utilities of Ann and Bob: a cloud of grey points for random allocations below a straight purple frontier; a dashed quarter circle shows the same frontier after relabelling utilities; the endowment lies inside and the market equilibrium lies on the frontier inside a short green segment](../../images/ch96-utility-frontier.svg)

*Figure 96.2 — The utility-possibility frontier. The market lands on it, but only somewhere in the short green segment the endowment allows. Its shape is not even meaningful: relabelling utilities (same preferences) bends a straight line into a circle.*

Figure 96.2 already hints at the hook's answer. Efficiency is a very weak demand. Every point of the frontier is efficient, including "Ann gets everything". Choosing *among* efficient points is a question of justice, and it requires comparing Ann's welfare with Bob's, a comparison that ordinal utility forbids.

### Arrow's theorem, proved

Let $A$ be a set of at least three alternatives and $N=\{1,\dots,n\}$ the voters. A **social welfare function** $F$ maps each profile of strict individual rankings to a strict social ranking. The axioms (stated in [§42.3](../c-game-theory/42-mechanism-design.md#423-the-mathematics)) are **unanimity** (if all rank $a$ above $b$, so does society) and **independence of irrelevant alternatives** (IIA: society's $a$-versus-$b$ verdict depends only on each voter's $a$-versus-$b$ verdict).

> [!theorem] Arrow (1951)
> If $\lvert A\rvert\ge3$ and $F$ satisfies unanimity and IIA, some voter $n^*$ is a dictator: society ranks $a$ above $c$ whenever $n^*$ does.

*Proof* (following Geanakoplos, 2005).

*Step 1: extremes stay extreme.* Suppose every voter puts alternative $b$ at the very top or very bottom of their ranking. We claim society does too. If not, there are $a,c$ with $a\succ b\succ c$ socially. Change every voter's ranking by moving $c$ just above $a$. No voter's $a$-versus-$b$ or $b$-versus-$c$ comparison changes, because $b$ is extreme. By IIA society still has $a\succ b$ and $b\succ c$, so $a\succ c$ by transitivity. But now every voter has $c$ above $a$, and unanimity demands $c\succ a$. Contradiction.

*Step 2: find the pivot.* Start with a profile where everyone ranks $b$ last, so society does (unanimity). Move $b$ to the top of voter $1$'s ranking, then voter $2$'s, and so on. At the end everyone ranks $b$ first and so does society. Let $n^*$ be the first voter whose move changes society's ranking of $b$. By Step 1 it jumps from bottom to top. Call the profiles just before and just after this move I and II.

*Step 3: $n^*$ dictates every pair not involving $b$.* Take $a,c\neq b$ and any profile in which $n^*$ ranks $a$ above $c$. Rearrange it, without changing anyone's $a$-versus-$c$ order, so that $n^*$ ranks $a\succ b\succ c$, voters before $n^*$ put $b$ first and voters after put $b$ last. On the pair $\{a,b\}$ every voter now agrees with profile I, where $b$ was socially last, so by IIA society has $a\succ b$. On $\{b,c\}$ every voter agrees with profile II, where $b$ was socially first, so $b\succ c$. Hence society has $a\succ c$. Everyone except $n^*$ ranked $a$ and $c$ arbitrarily, so by IIA society follows $n^*$ on $\{a,c\}$.

*Step 4: $n^*$ dictates pairs involving $b$ too.* Repeat Steps 2–3 with a third alternative $c$ in the role of $b$. This produces a voter who dictates every pair not involving $c$, including $\{a,b\}$. But going from profile I to profile II, only $n^*$ changed her ranking, and society's verdict on $a$ versus $b$ flipped. So nobody other than $n^*$ can dictate $\{a,b\}$, and the new dictator is $n^*$ again. ∎

**Kinship with Gibbard–Satterthwaite.** Given a strategy-proof voting rule $f$ that picks a winner and can elect anyone, define a social ranking by "$a$ above $b$ if $f$ elects $a$ when every voter lifts $a$ and $b$ (in their own order) to the top of their ballot". Strategy-proofness implies that this ranking satisfies IIA and is transitive, and ontoness gives unanimity. Arrow then supplies a dictator, and that dictator's favourite always wins under $f$. The pivotal voter of Step 2 is the "voter whose ranking flips society's choice" that [Ch. 42](../c-game-theory/42-mechanism-design.md) promised.

Why markets escape: an exchange economy isn't a social welfare function on arbitrary rankings. Preferences are restricted (more is better, convexity), the outcome is an allocation rather than a ranking, and prices use the *marginal rates* that IIA throws away. Restrict preferences for voting too, for example to **single-peaked** ones (each voter has an ideal point on a left–right line), and majority rule becomes transitive again (Black, 1948).

> [!mystery] How often can a voter cheat?
> Gibbard–Satterthwaite says every non-dictatorial rule is manipulable at *some* profile. Is manipulation common or vanishingly rare? Friedgut, Kalai and Nisan (2008), and then Mossel and Rácz (2015) for any number of alternatives, proved that for any rule far from every dictatorship, a random profile is manipulable with probability at least an inverse polynomial in $n$ and the number of alternatives. The proofs use Fourier analysis on the group $\{0,1\}^n$, which is the representation theory of [Ch. 95](95-representation-theory.md). What the *sharp* bound is remains open.

### Bargaining as a game: Rubinstein

Two players split a pie of size $1$. Player 1 offers a split, and player 2 accepts or rejects. After a rejection, one period passes and player 2 makes an offer, and so on forever. A payoff received $t$ periods from now is worth $\delta^t$ of its face value, with $0<\delta<1$.

> [!theorem] Rubinstein (1982)
> The game has a unique subgame-perfect equilibrium. The first proposer gets $\dfrac1{1+\delta}$, the responder accepts at once, and nothing is wasted.

*Proof* (Shaked and Sutton). Every subgame where a given player proposes is the same game. Let $M$ and $m$ be the supremum and infimum of the proposer's payoff over all subgame-perfect equilibria. A responder who rejects becomes the proposer next period, worth at most $\delta M$ today. So any offer giving her more than $\delta M$ is accepted, and the proposer can secure nearly $1-\delta M$: $m\ge1-\delta M$. Rejecting is worth at least $\delta m$ to the responder, so she never accepts less, and the proposer gets at most $1-\delta m$ (a rejected offer gives the proposer at most $\delta(1-m)\le1-\delta m$): $M\le1-\delta m$. Combining,

$$
M\le1-\delta m\le1-\delta(1-\delta M)\implies M(1-\delta^2)\le1-\delta\implies M\le\frac1{1+\delta},
$$

and then $m\ge1-\delta M\ge\frac1{1+\delta}$. So $M=m=\frac1{1+\delta}$. ∎

The first-mover advantage $\frac1{1+\delta}-\frac\delta{1+\delta}$ vanishes as offers come faster. With different patience rates, $\delta_i=e^{-r_i\Delta}$ for a period of length $\Delta$, the stationary equations $x_1=1-\delta_2x_2$, $x_2=1-\delta_1x_1$ give

$$
x_1=\frac{1-\delta_2}{1-\delta_1\delta_2}=\frac{1-e^{-r_2\Delta}}{1-e^{-(r_1+r_2)\Delta}}\;\longrightarrow\;\frac{r_2}{r_1+r_2}\quad(\Delta\to0).
$$

That is precisely the maximiser of the weighted Nash product $x^{1/r_1}(1-x)^{1/r_2}$. With equal patience it is the symmetric Nash solution $\tfrac12$ of [Ch. 44](../c-game-theory/44-coalitions-and-voting.md). Nash's axioms are what impatience looks like in the limit, and patience is bargaining power.

### Two more debts: Shapley and Myerson–Satterthwaite

**Why Shapley's axioms leave no freedom.** Every game $v$ on $N$ is uniquely a combination $v=\sum_{T\neq\varnothing}c_T\,u_T$ of **unanimity games** $u_T(S)=1$ if $T\subseteq S$, else $0$ (solve for the "dividends" $c_T$ by inclusion–exclusion). In $c\,u_T$, players outside $T$ are null and those inside are symmetric, so efficiency forces $c/\lvert T\rvert$ each. Additivity then fixes the value of every game: $\phi_i(v)=\sum_{T\ni i}c_T/\lvert T\rvert$. For the microgrid of Ch. 44 the dividends are $c_A=10$, $c_B=20$, $c_C=0$, $c_{AB}=10$, $c_{AC}=c_{BC}=40$, $c_{ABC}=-30$, and A gets $10+5+20-10=25$, as before. Drop additivity and you get the nucleolus, the Banzhaf index and other values, each defensible and each different.

**Why efficient trade can't balance its budget.** A seller with cost $c$ and a buyer with value $v$, both independent and uniform on $[0,1]$, want to trade whenever $v>c$. For any incentive-compatible mechanism that both agree to join, the payment identity of [Ch. 42](../c-game-theory/42-mechanism-design.md) (Myerson's lemma, applied to each side) shows that an exactly balanced budget requires expected **virtual surplus**

$$
E\Big[\big(\underbrace{2v-1}_{\text{buyer's virtual value}}-\underbrace{2c}_{\text{seller's virtual cost}}\big)\,q(v,c)\Big]\;\ge\;0,
$$

where $q$ is the probability of trade. For efficient trade, $q=\mathbf{1}\{v>c\}$, and the left side is $\int_0^1\int_0^v(2v-1-2c)\,dc\,dv=\int_0^1(v^2-v)\,dv=-\tfrac16<0$. So no such mechanism exists (Myerson and Satterthwaite, 1983). The best achievable rule trades only when $v-c\ge\tfrac14$. It makes the integral exactly $0$ and captures $\tfrac9{64}$ of the $\tfrac16$ available gains, about $84\%$.

## 96.4 The algorithm

Arrow's theorem is small enough to check by brute force. With two voters and three alternatives, IIA means a rule is just one function per pair from "(voter 1 prefers $x$?, voter 2 prefers $x$?)" to "society prefers $x$?". Unanimity fixes two of the four inputs, so each pair has two free bits and there are $2^6=64$ candidate rules. Keep those that never produce a cycle on any of the $6^2=36$ profiles:

```python
from itertools import permutations, product, combinations

alts = "abc"
rankings = list(permutations(alts))                  # 6 strict rankings
pairs = list(combinations(alts, 2))                  # (a,b), (a,c), (b,c)
prefers = lambda r, x, y: r.index(x) < r.index(y)

survivors = []
for bits in product([False, True], repeat=6):
    rule = {p: {(True, True): True, (False, False): False,
                (True, False): bits[2 * i], (False, True): bits[2 * i + 1]} for i, p in enumerate(pairs)}
    def social(profile, x, y):                       # does society rank x above y?
        if (x, y) in rule:
            return rule[(x, y)][tuple(prefers(r, x, y) for r in profile)]
        return not rule[(y, x)][tuple(prefers(r, y, x) for r in profile)]
    ok = all(not (social(P, x, y) and social(P, y, z) and social(P, z, x))
             for P in product(rankings, repeat=2) for x, y, z in permutations(alts))
    if ok:
        survivors.append(bits)
print(len(survivors), "transitive rules out of 64:", survivors)
```

Exactly two rules survive: "voter 1 decides every pair" and "voter 2 decides every pair". For three voters the same search has $2^{18}$ candidates, and it is still instant. Brute force grows doubly exponentially, though: $n$ voters give $2^{n}-2$ free bits per pair. Computer proofs of impossibility theorems (Tang and Lin, 2009, using a SAT solver) work by reducing to such a small base case and proving an induction step by hand.

Computing a market equilibrium is harder than checking one. The natural algorithm is **tâtonnement** ("groping"), Walras' own idea: raise the price of every good in excess demand, lower it for goods in excess supply, and repeat. Each step costs one evaluation of the excess-demand function $\mathbf{z}(\mathbf{p})=\sum_i\big(\mathbf{x}_i(\mathbf{p})-\boldsymbol\omega_i\big)$. For economies where goods are **gross substitutes** (raising one price raises demand for every other good), tâtonnement converges. In general it need not, and Scarf's algorithm (1967), which follows a path through a triangulated price simplex, was the first method guaranteed to find an approximate fixed point.

## 96.5 The experiment

First, drag the Edgeworth box. Ann's weight on apples is $0.6$ and Bob's is $0.4$; the endowment is $(8,2)$.

```desmos-graph
left=0; right=10; bottom=0; top=10
---
x^{0.6}y^{0.4}=8^{0.6}\cdot2^{0.4}|blue
(10-x)^{0.4}(10-y)^{0.6}=2^{0.4}\cdot8^{0.6}|red
0.6(1-0.4)y(10-x)=0.4(1-0.6)x(10-y)|purple
x+y=10|black|dashed
(8,2)|label:endowment
(6,4)|label:equilibrium
```

*Ann's (blue) and Bob's (red) indifference curves through the endowment, the contract curve (purple) and the budget line at price ratio 1.*

Change the endowment $(8,2)$ to $(9,5)$ (replace both $8^{0.6}\cdot2^{0.4}$ and $2^{0.4}\cdot8^{0.6}$ by the new holdings, Ann's in the first line and Bob's $(1,5)$ in the second). The lens moves, but the contract curve does not, because it depends only on preferences. Then change the weights $0.6$ and $0.4$ to equal values $0.5$ and watch the contract curve straighten into the diagonal.

Then run tâtonnement, first on the box economy and then on Herbert Scarf's 1960 economy. In Scarf's economy three consumers each own one good and want two goods in fixed proportions.

```python
import numpy as np

# Tatonnement: raise the price of whatever is over-demanded, lower it if under-demanded.
a, b = 0.6, 0.4                 # Cobb-Douglas weights on apples for Ann and Bob
wA = np.array([8.0, 2.0])       # Ann's endowment (apples, bread); Bob has the rest
wB = np.array([10.0, 10.0]) - wA

def excess_box(p):              # p = (p_apple, p_bread)
    IA, IB = p @ wA, p @ wB
    xA = np.array([a * IA / p[0], (1 - a) * IA / p[1]])
    xB = np.array([b * IB / p[0], (1 - b) * IB / p[1]])
    return xA + xB - wA - wB

def excess_scarf(p):            # Scarf (1960): consumer i owns good i, wants goods i and i+1 in fixed 1:1 ratio
    z = -np.ones(3)
    for i in range(3):
        k = (i + 1) % 3
        q = p[i] / (p[i] + p[k])          # consumer i's income buys q units of each of goods i and k
        z[i] += q; z[k] += q
    return z

def tatonnement(excess, p, h=0.05, steps=4000):
    path = [p.copy()]
    for _ in range(steps):
        p = np.maximum(p + h * excess(p), 1e-9)
        p = p / p.sum()                   # only relative prices matter
        path.append(p.copy())
    return np.array(path)

box = tatonnement(excess_box, np.array([0.8, 0.2]))
print("box economy: final prices", box[-1].round(4), " excess demand", excess_box(box[-1]).round(6))

scarf = tatonnement(excess_scarf, np.array([0.5, 0.3, 0.2]))
dist = np.linalg.norm(scarf - 1 / 3, axis=1)
print("Scarf economy: distance from equilibrium (1/3,1/3,1/3) at steps 0, 1000, 2000, 4000:",
      dist[[0, 1000, 2000, 4000]].round(4))
```

The box economy settles at equal prices $(0.5,0.5)$, our $p=1$. Scarf's economy has a unique equilibrium at equal prices, yet the distance from it *grows* from $0.22$ to about $0.5$, and prices circle the equilibrium without ever approaching it. Things to try: shrink the step `h` to $0.005$ (does the spiral stop growing?), start Scarf's prices at exactly $(1/3,1/3,1/3)$ plus a tiny nudge, and plot `scarf[:, 0]` against `scarf[:, 1]` to see the orbit.

> [!mystery] The market that never settles
> In continuous time Scarf's prices move on closed orbits around the equilibrium forever. His economy is perfectly well behaved (unique equilibrium, continuous demand), yet the invisible hand never finds it. Whether an equilibrium is *stable* is a question about a dynamical system $\dot{\mathbf{p}}=\mathbf{z}(\mathbf{p})$, and centres, spirals and limit cycles are classified in [Ch. 83](../i-dynamical-systems/83-dynamical-systems.md).

## 96.6 The challenge

1. **Intuition.** In the Edgeworth box, why must the contract curve pass through both corners? Is the equilibrium's position on the contract curve determined by preferences, by endowments, or by both?
2. **Derivation.** Borda count (Ch. 44) satisfies unanimity and isn't a dictatorship, so by Arrow it violates IIA. Construct an explicit three-voter, three-candidate pair of profiles in which no voter changes her $a$-versus-$b$ order but Borda's social $a$-versus-$b$ order flips.
3. **Implementation.** Extend the brute-force Arrow search to three voters ($2^{18}$ rules). Then drop unanimity (all $4$ inputs free per pair, $2^{24}$ rules for three voters) and list the transitive IIA rules that survive. Wilson's theorem (1972) predicts what they are: describe them.
4. **Open investigation.** Run tâtonnement on random exchange economies with three goods and Cobb–Douglas, Leontief or CES consumers. Estimate how often it fails to converge as preferences move from substitutes towards complements. No general characterisation of the economies on which simple price dynamics converge is known.

## 96.7 Loose threads

The hook asked how markets can aggregate preferences perfectly while votes cannot aggregate them at all. Now we can say. The first welfare theorem asks for very little (efficiency only, which even "Ann gets everything" satisfies), and it gets to use prices, which encode marginal rates of substitution. Arrow asks for much more, a full social ranking that works on every possible profile, while IIA forbids it to use intensity. Prices carry information that ballots throw away, and markets promise something that ballots are asked to exceed. Rubinstein added that even raw haggling, given impatience, ends at Nash's axiomatic split. Myerson–Satterthwaite added a limit on what any market can do when people hold private information.

What we didn't tell you:

- **Common values and the winner's curse.** In the auctions of [Ch. 43](../c-game-theory/43-auctions-and-markets.md), when an item is worth the same unknown amount to everyone, winning means you were the most optimistic. Milgrom and Weber (1982) made "affiliated" signals precise: high signals make other high signals more likely. They proved the **linkage principle**, which says that formats tying the price more closely to other bidders' information raise expected revenue. So the English auction beats the second-price auction, which beats the first-price. The 2020 Nobel prize was awarded for the spectrum auctions built on this.
- **Which equilibrium is reasonable?** The Intuitive Criterion of [Ch. 41](../c-game-theory/41-incomplete-information-games.md) kills implausible beliefs in signalling games, but stronger refinements (divinity, stability in the sense of Kohlberg and Mertens) can disagree with it, and experiments find people often play the "unreasonable" pooling equilibria. Whether refinement theory describes behaviour or only consistency is still argued.
- **Large markets.** Debreu and Scarf (1963) proved that as an economy is replicated, the **core** (allocations no coalition can block) shrinks to the set of Walrasian equilibria. With many traders, competitive prices are the only allocations that survive coalitional bargaining.

> [!mystery] Anything can happen, and finding it is hard
> Individual demand functions have strong structure, but the Sonnenschein–Mantel–Debreu theorem (1972–74) says aggregate excess demand can be essentially *any* continuous function satisfying Walras' law. So economic theory alone predicts almost nothing about the number, location or stability of equilibria. And finding one is, in general, as hard as finding a Brouwer fixed point: computing equilibria of exchange economies is PPAD-complete. What PPAD is, and why guaranteed existence doesn't make search easy, is in [Ch. 30](../a-algorithms/30-limits-of-computation.md).

Each of these is a doorway into the economics research literature. The proofs you have just seen, an inequality summed over consumers, a pivot hunted down by moving one voter at a time, a sup and an inf squeezed together, are typical of how that literature works: short arguments around a clever choice of what to compare. How such arguments are found is the subject of [Ch. 98](98-mathematical-research.md).

## Further reading

- A. Mas-Colell, M. D. Whinston, J. R. Green, *Microeconomic Theory*, Oxford University Press, 1995. Chapters 1–3, 15–17 and 21–23 cover everything here at graduate level.
- K. J. Arrow, *Social Choice and Individual Values*, 3rd ed., Yale University Press, 2012. The original thesis, with later reflections.
- J. Geanakoplos, "Three brief proofs of Arrow's impossibility theorem", *Economic Theory* 26 (2005) 211–215. The proof given above, and two others.
- M. J. Osborne and A. Rubinstein, *Bargaining and Markets*, Academic Press, 1990 (free online). Rubinstein's model and its connection to Nash and to markets.
- H. Moulin, *Fair Division and Collective Welfare*, MIT Press, 2003. Shapley values, bargaining and fair division in one framework.
