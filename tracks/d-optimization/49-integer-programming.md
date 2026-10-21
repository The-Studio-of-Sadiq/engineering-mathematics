# Chapter 49 — Integer Programming and Combinatorial Optimization

[← Previous](48-duality.md) · [Track D](README.md) · [Table of Contents](../../00-front/table-of-contents.md) · [Next →](50-operations-research.md)

> **Level 2 · Track D — Optimization and Operations Research**
> **Builds on:** [Ch. 8](../../part-3-linear-algebra/08-matrices-and-linear-systems.md), [Ch. 23](../../part-7-numerical-computational/23-optimization.md), [Ch. 47](47-convexity.md), [Ch. 48](48-duality.md) · **Leads to:** [Ch. 50](50-operations-research.md), [Ch. 52](52-dynamic-programming.md) · **Related:** [Ch. 27](../a-algorithms/27-graph-theory.md), [Ch. 28](../a-algorithms/28-combinatorics.md), [Ch. 30](../a-algorithms/30-limits-of-computation.md)
> **Reading time:** ~30 min

## 49.1 The hook

A data-centre team is filling a hall with server racks. There are two designs. A type-X rack delivers 5 units of compute and draws 5 kW; a type-Y rack delivers 8 units and draws 9 kW. The hall has 6 rack positions and a 45 kW power feed. How many of each should they install?

An engineer who knows linear programming ([§23.3](../../part-7-numerical-computational/23-optimization.md#233-linear-programming)) writes it down in a minute:

$$
\max\ 5x+8y\quad\text{s.t.}\quad x+y\le6,\quad 5x+9y\le45,\quad x,y\ge0,
$$

and the solver replies at once: $x=2.25$, $y=3.75$, compute $41.25$.

You can't install a quarter of a rack. The obvious fix is to round. Rounding to the nearest integers gives $(2,4)$, which draws $10+36=46$ kW and trips the breaker. Rounding down to $(2,3)$ is safe but gives only 34, 18% below the LP's promise. Hunting around nearby, $(3,3)$ gives 39. The true optimum turns out to be $(0,5)$: five type-Y racks, compute 40, and *zero* of the type the LP wanted 2.25 of. The best integer answer isn't a rounding of the LP answer at all. It sits on the far side of the feasible region.

With two variables you can list all 25 feasible integer points and check. Real problems aren't like that. An airline assigning crews to flights, a chip designer placing logic cells, a utility deciding which power stations to switch on tomorrow: each has thousands or millions of yes/no decisions. One hundred binary choices already give $2^{100}\approx1.27\times10^{30}$ combinations. At a billion checks per second that would take about $4\times10^{13}$ years, thousands of times the age of the universe.

Yet commercial solvers settle problems like these to proven optimality every day, often in minutes. They don't round, and they don't enumerate. This chapter is about the trick: **use the easy continuous problem as a bound, and cut the hopeless search tree down with it.**

## 49.2 The big idea

The integer points sit inside the LP's polygon like seeds in a fruit. The LP relaxation, which drops the integrality requirement, is a *convex* problem ([Ch. 47](47-convexity.md)): easy to solve, and it always over-promises. Its value is an upper bound on what any integer design can achieve. Integer programming runs on two ideas built on that bound.

**Idea 1: divide and bound.** Split the problem into pieces, "$x\le2$" and "$x\ge3$", which together contain every integer point and remove the fractional LP solution in between. Solve each piece's LP. If a piece's bound is already worse than a design you hold in hand, the whole piece, with every integer point inside it, can be thrown away without being looked at. That is **branch-and-bound**. Its power comes entirely from how much the bounds let you skip.

**Idea 2: tighten the fruit.** Between the integer points and the LP polygon there's wasted space. If you could shrink the polygon to the **convex hull** of the integer points, the tightest convex shape containing them, its corners would *be* integer points, and a single LP would solve the integer problem. You can't compute the hull directly; it may have astronomically many facets. But you can add a few of its walls, **cutting planes**, exactly where the LP solution is cheating.

![The LP polygon with integer lattice points; the LP optimum (2.25,3.75) is a fractional corner; rounding it gives an infeasible point (2,4) or a poor point (2,3); the integer optimum is (0,5); the convex hull of integer points is shown dashed, with the cut 2x+3y<=15 as one of its edges](../../images/ch49-relaxation.svg)

*The rack problem. The LP relaxation (blue) has a fractional corner worth 41.25. The convex hull of the integer points (green dashed) is smaller. One extra wall, the cut $2x+3y\le15$ (purple), slices off the fractional corner, and on the trimmed polygon the LP optimum is the integer point $(0,5)$.*

A useful analogy is detective work. Branching is splitting the suspects into groups and interrogating each group's alibi. A bound is an alibi so strong that the whole group goes home at once. A cut is a new fact that rules out an impossible story ("nobody was in two rooms at the same time") without naming a suspect. Real solvers do both, constantly. That is **branch-and-cut**.

## 49.3 The mathematics

### Modelling with integers

The step from "continuous" to "integer" is mostly about what you can *say*. Binary variables $z\in\{0,1\}$ turn logic into linear inequalities.

| Requirement | Constraint |
|---|---|
| choose at most one of items 1–4 | $z_1+z_2+z_3+z_4\le1$ |
| if you build A you must build B | $z_A\le z_B$ |
| power $p$ only if the station is on, $0\le p\le P_{\max}$ | $0\le p\le P_{\max}z$ |
| fixed start-up cost $F$ plus running cost $c$ | cost $Fz+cp$ with $p\le P_{\max}z$ |
| either $\mathbf{a}^\mathsf{T}\mathbf{x}\le b_1$ or $\mathbf{c}^\mathsf{T}\mathbf{x}\le b_2$ | $\mathbf{a}^\mathsf{T}\mathbf{x}\le b_1+Mz$, $\mathbf{c}^\mathsf{T}\mathbf{x}\le b_2+M(1-z)$ |

The last row uses a **big-$M$**, a constant large enough that the constraint switches off when its binary says so. The fixed-charge row is exactly the start-up cost that broke the electricity market's prices in [§48.3](48-duality.md#483-the-mathematics). Non-convexity that defeats prices can often be *modelled* exactly with a few binaries.

The classic test problem is the **0–1 knapsack**: $\max\sum v_iz_i$ s.t. $\sum w_iz_i\le W$, $z_i\in\{0,1\}$. Pack items of value $v_i$ and weight $w_i$ into a bag of capacity $W$. It looks innocent. It is NP-hard.

> [!warning] Common Pitfall: big-$M$ too big
> Any $M$ large enough gives a *correct* model, but not an equally good one. In the relaxation, $z$ may take the value $p/P_{\max}$, so a huge $M$ lets $z$ be tiny and almost cancels the fixed cost. The LP bound becomes weak, and branch-and-bound must work far harder. Choose the smallest $M$ that is valid. Physical limits ($P_{\max}$ itself) are usually the right choice.

### Relaxations bound, and hulls are exact

Write the integer program as $z_{\mathrm{IP}}=\max\{\mathbf{c}^\mathsf{T}\mathbf{x}:\mathbf{x}\in S\}$ with $S=P\cap\mathbb{Z}^n$, where $P=\{\mathbf{x}:\mathbf{A}\mathbf{x}\le\mathbf{b}\}$ is the LP polyhedron, and let $z_{\mathrm{LP}}=\max\{\mathbf{c}^\mathsf{T}\mathbf{x}:\mathbf{x}\in P\}$.

> [!theorem] Theorem 49.1 — The relaxation bound and the convex hull
> 1. $z_{\mathrm{IP}}\le z_{\mathrm{LP}}$.
> 2. If $S$ is finite and nonempty, $z_{\mathrm{IP}}=\max\{\mathbf{c}^\mathsf{T}\mathbf{x}:\mathbf{x}\in\operatorname{conv}S\}$, and the maximum is attained at a point of $S$.

*Proof.* (1) $S\subseteq P$, and maximising over a larger set can't give less. (2) Every $\mathbf{x}\in\operatorname{conv}S$ is a convex combination $\sum\theta_k\mathbf{s}_k$ of points of $S$ ([§47.3](47-convexity.md#473-the-mathematics)), so $\mathbf{c}^\mathsf{T}\mathbf{x}=\sum\theta_k\mathbf{c}^\mathsf{T}\mathbf{s}_k\le\max_k\mathbf{c}^\mathsf{T}\mathbf{s}_k\le z_{\mathrm{IP}}$. A linear function can't do better on a mixture than on the best ingredient. The reverse inequality holds because $S\subseteq\operatorname{conv}S$. $\square$

So integer programming *is* linear programming over the right polygon, $\operatorname{conv}S$. The difficulty is entirely in describing that polygon. For the racks, $\operatorname{conv}S$ has corners $(0,0)$, $(6,0)$, $(3,3)$, $(0,5)$ and the inequality description

$$
x\ge0,\qquad y\ge0,\qquad x+y\le6,\qquad 2x+3y\le15 .
$$

The **integrality gap** $z_{\mathrm{LP}}-z_{\mathrm{IP}}=41.25-40=1.25$ measures how much the relaxation over-promises. It's the integer-programming cousin of the duality gap in [Ch. 48](48-duality.md).

### Branch-and-bound

Branch-and-bound explores a tree. Each node is the original problem plus some extra bounds on variables, and its LP value $\bar z$ is an upper bound for every integer point in that node. The best integer solution found so far, the **incumbent**, has value $z^\circ$.

> [!theorem] Theorem 49.2 — Pruning never loses the optimum
> Suppose the integer points of a node are split into two children by $x_j\le\lfloor\bar x_j\rfloor$ and $x_j\ge\lfloor\bar x_j\rfloor+1$, and a node is discarded only if (a) its LP is infeasible, (b) its LP solution is integer (it then updates the incumbent), or (c) $\bar z\le z^\circ$. Then when no open nodes remain, the incumbent is optimal.

*Proof.* Every integer point of a node lies in exactly one child, because no integer lies strictly between $\lfloor\bar x_j\rfloor$ and $\lfloor\bar x_j\rfloor+1$. So the leaves of the tree always partition $S$. Take an optimal $\mathbf{x}^\star$ and follow it down the tree to the leaf that contains it. That leaf was closed for one of three reasons. (a) is impossible, since $\mathbf{x}^\star$ is feasible for it. Under (b), the leaf's LP optimum is an integer point of value at least $\mathbf{c}^\mathsf{T}\mathbf{x}^\star$, and it entered the incumbent. Under (c), $\mathbf{c}^\mathsf{T}\mathbf{x}^\star\le\bar z\le z^\circ$. In each case the final incumbent is at least as good as $\mathbf{x}^\star$. The tree is finite when the variables are bounded, because each branch shrinks an integer range. $\square$

Here is the whole search for the racks, branching on the more fractional variable and exploring the "$\le$" child first.

![The branch-and-bound tree: root LP (2.25,3.75) value 41.25 splits on x into x<=2 with LP 41.11 and x>=3 with LP (3,3) value 39; the left node splits on y into an integer leaf (2,3) value 34 and a node (1.8,4) value 41, which splits on x into an infeasible leaf and a node (1,4.44) value 40.56, which splits on y into an integer leaf (1,4) value 37 and the integer optimum (0,5) value 40](../../images/ch49-branch-and-bound.svg)

*Branch-and-bound on the rack problem. Each box shows the node's LP solution and its value, an upper bound for that branch. Depth-first search finds the incumbent 34, improves it to 40 at $(0,5)$, and then closes the right branch without branching on it: that branch's integer optimum is only 39.*

Nine LPs settle the problem, and the final one is closed by its bound. In larger problems the bound does far more work: whole subtrees with billions of leaves vanish in a single comparison.

### Cutting planes

> [!theorem] Theorem 49.3 — Chvátal–Gomory cuts are valid
> Let $\mathbf{A}\mathbf{x}\le\mathbf{b}$, $\mathbf{x}\ge\mathbf{0}$ define $P$ and let $\mathbf{u}\ge\mathbf{0}$. Then every integer point of $P$ satisfies
> $$
> \sum_j\big\lfloor(\mathbf{u}^\mathsf{T}\mathbf{A})_j\big\rfloor x_j\ \le\ \big\lfloor\mathbf{u}^\mathsf{T}\mathbf{b}\big\rfloor .
> $$

*Proof.* Combining the rows with weights $\mathbf{u}\ge\mathbf{0}$ gives the valid inequality $\sum_j(\mathbf{u}^\mathsf{T}\mathbf{A})_jx_j\le\mathbf{u}^\mathsf{T}\mathbf{b}$, the same "combine constraints" move as in Farkas' lemma and LP duality. Because $\mathbf{x}\ge\mathbf{0}$, rounding the coefficients *down* can only shrink the left side, so $\sum_j\lfloor(\mathbf{u}^\mathsf{T}\mathbf{A})_j\rfloor x_j\le\mathbf{u}^\mathsf{T}\mathbf{b}$ still holds. For integer $\mathbf{x}$ the left side is now an integer, and an integer at most $\mathbf{u}^\mathsf{T}\mathbf{b}$ is at most $\lfloor\mathbf{u}^\mathsf{T}\mathbf{b}\rfloor$. $\square$

The last step is the only place integrality enters, and it's where fractional points get cut off. For the racks take $\mathbf{u}=(\tfrac34,\tfrac14)$:

$$
\tfrac34(x+y)+\tfrac14(5x+9y)\le\tfrac34\cdot6+\tfrac14\cdot45\quad\Longrightarrow\quad 2x+3y\le15.75\quad\Longrightarrow\quad 2x+3y\le15 .
$$

The LP optimum has $2(2.25)+3(3.75)=15.75$, so it violates the cut by $0.75$. Re-solving with the cut added gives $(0,5)$, value 40, the integer optimum, in a single round. The cut is exactly the missing facet of the hull.

That was lucky. In general Gomory (1958) showed how to read a violated cut directly from the final simplex tableau. Chvátal (1973) proved that repeating the round-down operation finitely many times always reaches $\operatorname{conv}S$, but the number of rounds can be enormous. In practice cuts and branching work together: cut until progress stalls, then branch.

### When the LP is already integral

Some integer programs need no integer machinery at all.

> [!definition] Total unimodularity
> A matrix $\mathbf{A}$ is **totally unimodular** (TU) if every square submatrix has determinant $0$, $+1$ or $-1$.

> [!theorem] Theorem 49.4 — TU matrices give integer vertices
> If $\mathbf{A}$ is TU and $\mathbf{b}$ is an integer vector, every vertex of $P=\{\mathbf{x}:\mathbf{A}\mathbf{x}\le\mathbf{b},\ \mathbf{x}\ge\mathbf{0}\}$ is integer. Hence the LP relaxation solves the integer program.

*Proof.* Stacking $\mathbf{A}$ on $-\mathbf{I}$ (for $\mathbf{x}\ge\mathbf{0}$) keeps the matrix TU, because expanding a determinant along the rows coming from $-\mathbf{I}$ leaves a $\pm$ determinant of a submatrix of $\mathbf{A}$. A vertex is the unique solution of $\mathbf{B}\mathbf{x}=\mathbf{b}'$ for some nonsingular square submatrix $\mathbf{B}$ of the stacked matrix and the matching integer entries $\mathbf{b}'$. By Cramer's rule ([Ch. 8](../../part-3-linear-algebra/08-matrices-and-linear-systems.md)), $x_i=\det\mathbf{B}_i/\det\mathbf{B}$, where $\mathbf{B}_i$ is $\mathbf{B}$ with column $i$ replaced by $\mathbf{b}'$. The numerator is an integer, and the denominator is $\pm1$. $\square$

Incidence matrices of directed graphs, and of bipartite graphs, are TU. That's why shortest paths, maximum flows and the assignment problem ($n$ workers, $n$ jobs) are "integer programs" that a plain LP, or a specialised network algorithm ([Ch. 27](../a-algorithms/27-graph-theory.md)), solves exactly in polynomial time. Recognising hidden network structure in an integer model is one of the most profitable skills in operations research.

### Lagrangian relaxation

Duality ([Ch. 48](48-duality.md)) gives a second kind of bound. Move the awkward constraints into the objective with prices $\boldsymbol\lambda\ge\mathbf{0}$ and keep the easy ones, integrality included. The dual function is then an upper bound on $z_{\mathrm{IP}}$ (for a maximisation), and a classical result of Geoffrion (1974) says the best such bound is *never worse* than the LP bound, and often strictly better. The gap of Worked Example 48.1 doesn't go away, but it can be pushed below the integrality gap. For the crew-scheduling and unit-commitment problems of [Ch. 50](50-operations-research.md), that is often the difference between a bound you can use and one you can't.

### How hard is this, really?

Deciding whether a 0–1 knapsack can reach value $V$ is NP-complete, and general integer programming is NP-hard. Unless P = NP, no algorithm solves every instance in time polynomial in the input size. Branch-and-bound offers no guarantee either: there are tiny-looking families on which it must explore exponentially many nodes (§49.4 and §49.5 below show one).

> [!mystery] Hard in theory, easy in practice?
> Solvers routinely close industrial integer programs with millions of variables, yet the theory says the worst case is hopeless. Is that because real instances have hidden structure, or because P = NP after all? Nobody knows. The P versus NP question is open, and [Ch. 30](../a-algorithms/30-limits-of-computation.md) explains what a proof either way would have to overcome. Even the practical success of branch-and-cut has no satisfying theoretical explanation.

## 49.4 The algorithm

### Branch-and-bound, as code

```text
input: max c.x  s.t.  A x <= b,  x >= 0,  x_j integer for j in J
incumbent = None, z_best = -infinity
open = [ root problem ]
while open is not empty:
    node = select(open)                    # depth-first, best-bound, or a mix
    solve the LP relaxation of node -> (x_bar, z_bar)
    if infeasible or z_bar <= z_best:  continue          # prune
    if x_bar_j is integer for all j in J:
        incumbent, z_best = x_bar, z_bar;  continue      # new incumbent
    j = choose a fractional variable       # branching rule
    open += [ node + {x_j <= floor(x_bar_j)},  node + {x_j >= floor(x_bar_j) + 1} ]
report incumbent, and the gap  (max bound over open nodes) - z_best
```

Three choices decide whether this finishes in a second or in a century.

- **Node selection.** *Depth-first* dives quickly to integer leaves, so it finds incumbents early and uses little memory. *Best-bound* always expands the node with the highest $\bar z$, so it proves optimality with the fewest nodes but can keep a huge open list. Modern solvers dive first and switch later.
- **Branching variable.** "Most fractional" is the textbook rule and is often poor. *Strong branching* tentatively solves both children's LPs for several candidates and picks the variable whose children improve the bound most. It's expensive per node, but it often cuts the tree by orders of magnitude.
- **Warm starts.** A child differs from its parent by one bound, so the parent's optimal basis is still dual-feasible. A few **dual simplex** pivots ([§23.3](../../part-7-numerical-computational/23-optimization.md#233-linear-programming)) re-solve it, far cheaper than starting from scratch.

**Anytime answers.** At every moment the solver holds an incumbent (a lower bound) and the best open-node bound (an upper bound). The **optimality gap** $(\bar z_{\max}-z^\circ)/\lvert z^\circ\rvert$ tells you how far from optimal the current design could be. Industrial users often stop at 1%: a guarantee that rounding never gave them.

**Branch-and-cut** inserts a cut-separation loop at each node: solve the LP, look for violated cuts (Gomory, knapsack covers, clique inequalities…), add them, re-solve, and branch only when the cuts stop helping. Add primal heuristics that try to round or repair LP solutions into good incumbents, and you have the skeleton of every commercial MIP solver.

**Failure cases.** The worst case is exponential, and it isn't exotic. Take $2z_1+2z_2+\dots+2z_n=n$ with $n$ odd and $z_i\in\{0,1\}$. No solution exists (the left side is even), but setting every free variable to $\tfrac12$-ish keeps the LP feasible until about half the variables are fixed. So plain branch-and-bound must explore at least about $2^{n/2}$ nodes before proving infeasibility. A single cut ($\sum z_i\le\lfloor n/2\rfloor$ from rounding) or a parity argument ends it at once, which shows why cuts matter. Other trouble spots are weak big-$M$ formulations (§49.3), symmetry (identical machines make many branches equivalent copies of each other) and numerical tolerance ("integer" means within $10^{-6}$ or so, which can mislead with huge coefficients).

## 49.5 The experiment

**Slide the objective.** The graph draws the rack constraints and the integer points that matter. The dashed objective line $5x+8y=41.25$ touches the LP corner. Change `41.25` to `40`: the line drops to the integer optimum $(0,5)$, and no integer point lies between the two lines. Now change the objective to $5x+7y$ (edit both `8`s and recompute the level). Does the integer optimum stay at $(0,5)$?

```desmos-graph
left=-0.5; right=7; bottom=-0.5; top=6.5
---
y=6-x|blue|0<x<6
y=\frac{45-5x}{9}|blue|0<x<6
y=\frac{15-2x}{3}|purple|dotted|0<x<7.5
y=\frac{41.25-5x}{8}|red|dashed
(2.25,3.75)|label:LP optimum
(0,5)|label:integer optimum
(3,3)
(2,3)
(2,4)
```

*The rack problem in Desmos. Blue: the LP constraints. Purple dotted: the Chvátal–Gomory cut. Red dashed: an objective level line to slide.*

**Watch the tree explode.** The script solves a 0–1 knapsack by depth-first branch-and-bound. The bound at each node is the knapsack's LP relaxation, which has a closed form: sort items by value per weight, take them greedily, and take a fraction of the first item that doesn't fit (Dantzig's bound). The script needs only numpy.

```python
import numpy as np
import matplotlib.pyplot as plt

# --- knobs -------------------------------------------------------------
n, seed = 30, 1                       # number of items, random instance
rng = np.random.default_rng(seed)
w = rng.integers(10, 100, n)          # weights   (try: v = w + 10 for a hard instance)
v = rng.integers(10, 100, n)          # values
cap = int(0.5 * w.sum())              # capacity

order = np.argsort(-v / w); v, w = v[order], w[order]   # best value per weight first

def lp_bound(k, room, val):           # LP relaxation of items k..n-1: greedy + one fraction
    for i in range(k, n):
        if w[i] <= room: room -= w[i]; val += v[i]
        else: return val + v[i] * room / w[i]
    return val

best, nodes, trace = 0, 0, []
stack = [(0, cap, 0)]                 # (next item, remaining room, value so far)
while stack:
    k, room, val = stack.pop(); nodes += 1
    best = max(best, val)
    if k == n or lp_bound(k, room, val) <= best:
        continue                      # leaf, or pruned: this branch cannot beat 'best'
    stack.append((k + 1, room, val))                     # skip item k
    if w[k] <= room:
        stack.append((k + 1, room - w[k], val + v[k]))   # take item k (explored first)
    trace.append(best)

print(f"optimum {best},  root LP bound {lp_bound(0, cap, 0):.2f},  nodes {nodes}"
      f"  (full tree would have {2**(n+1)-1:.3g})")
plt.plot(trace); plt.xlabel("nodes branched"); plt.ylabel("best value found")
plt.title("incumbent over the search"); plt.show()
```

With the default random instance ($n=30$) the optimum is $1275$ (checked against an exact dynamic program), the root bound is $1282.85$, and the search visits only **59 nodes** out of a full tree of about $2.15\times10^9$. The bound is so tight that almost everything is pruned. Doubling to $n=50$ still needs only 169 nodes.

Now make the instance nasty. Replace the three instance lines with `w = 2 * rng.integers(5, 50, n)`, `v = w.copy()` and `cap = int(0.5 * w.sum()) | 1`. Every weight is even and the capacity is odd, so the bag can never be filled exactly, and the LP bound always sits one unit above anything achievable. Nothing gets pruned early. The node counts become:

| $n$ | 16 | 20 | 24 |
|---|---|---|---|
| nodes visited | 28,543 | 222,024 | 8,475,918 |

That's roughly an eight- to forty-fold increase for every four extra items: exponential growth, live on your laptop. (The $n=24$ run can take a few minutes.) Things to try: switch the stack to a best-bound priority queue and compare node counts; add the parity cut "total weight $\le$ cap $-1$" to `lp_bound` and watch the nasty instance collapse.

## 49.6 The challenge

1. **Intuition.** Draw a two-variable integer program whose LP relaxation is feasible and bounded but which has **no** integer feasible point at all. Then draw one where the LP optimum is far from *every* integer optimum, at least 10 units away.
2. **Derivation.** For the rack problem, derive the cut $2x+3y\le15$ as a Gomory fractional cut from the optimal simplex tableau of the LP relaxation (introduce slacks, find the row of a fractional basic variable, and take fractional parts). Check that you obtain the same inequality as the Chvátal–Gomory weights $\mathbf{u}=(\tfrac34,\tfrac14)$.
3. **Implementation.** Write a general branch-and-bound for small integer programs in $\max\{\mathbf{c}^\mathsf{T}\mathbf{x}:\mathbf{A}\mathbf{x}\le\mathbf{b},\ \mathbf{x}\in\mathbb{Z}^n_{\ge0}\}$ form, using `scipy.optimize.linprog` for the node LPs. Reproduce the tree of Figure 49.2. Then compare depth-first and best-bound node counts on random instances with $n=10$ to $20$.
4. **Open investigation.** For random knapsack instances (values and weights uniform on $[1,R]$), how does the expected number of branch-and-bound nodes grow with $n$? Polynomially? Does the answer change if values and weights are correlated, as in the `v = w + 10` variant? Plot node counts on log–log and log–linear axes and argue from your data which growth law fits.

## 49.7 Loose threads

Integer programming is where convexity ends and cleverness begins. Three of the deepest questions it raises are answered elsewhere.

> [!mystery] A different way to slice the knapsack
> Branch-and-bound split the knapsack by *items*. There's another way to organise the same search, by *remaining capacity*: "what is the best value achievable with capacity $c$ using the first $k$ items?" That table has only $nW$ entries, so it solves the knapsack in time proportional to $n$ times the *number* $W$, which can be fast or hopeless depending on how big $W$ is. Why doesn't this contradict NP-hardness, and what principle makes the table correct? Bellman's principle of optimality in [Ch. 52](52-dynamic-programming.md).

> [!mystery] Models too big to write down
> An airline's crew-scheduling model has one binary variable for every *legal pairing*, a sequence of flights a crew can fly in one trip. There are billions of them, far too many to list, let alone branch on. Yet such models are solved routinely by generating variables only when the dual prices of [Ch. 48](48-duality.md) say a new one would help. How does an LP solver "price out" a column it has never seen? That is column generation, in [Ch. 50](50-operations-research.md).

The third thread is the P versus NP mystery of §49.3 ([Ch. 30](../a-algorithms/30-limits-of-computation.md)). For counting how many integer points a polytope holds, and why that number grows like a polynomial when the polytope is scaled up, see the combinatorics of [Ch. 28](../a-algorithms/28-combinatorics.md).

## Further reading

- L. A. Wolsey, *Integer Programming*, 2nd ed., Wiley, 2020. The standard introduction: formulations, relaxations, branch-and-bound and cutting planes.
- M. Conforti, G. Cornuéjols and G. Zambelli, *Integer Programming*, Springer, 2014. Rigorous polyhedral theory, including total unimodularity and the Chvátal closure.
- A. Schrijver, *Theory of Linear and Integer Programming*, Wiley, 1986. The reference for the deep theory.
- F. S. Hillier and G. J. Lieberman, *Introduction to Operations Research*, McGraw-Hill (many editions). Accessible modelling examples, including the classic "rounding fails" illustration.
