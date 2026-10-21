import sys, pathlib; sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
from style import plt, np, PALETTE, INK, MUTED, new_figure, diagram, node, box, arrow, save
from itertools import permutations

# ---------------------------------------------------------------- Shapley value bars
v = {frozenset(): 0, frozenset("A"): 10, frozenset("B"): 20, frozenset("C"): 0,
     frozenset("AB"): 40, frozenset("AC"): 50, frozenset("BC"): 60, frozenset("ABC"): 90}
players = "ABC"
phi = {p: 0.0 for p in players}
for order in permutations(players):
    s = frozenset()
    for p in order:
        phi[p] += (v[s | {p}] - v[s]) / 6
        s = s | {p}
alone = [v[frozenset(p)] for p in players]
schemes = {
    "stand-alone value": alone,
    "equal split": [30, 30, 30],
    "proportional to stand-alone": [90 * a / sum(alone) for a in alone],
    "Shapley value": [phi[p] for p in players],
}
fig, ax = new_figure(figsize=(7.4, 4.3))
x = np.arange(3)
w = 0.2
cols = [MUTED, PALETTE[5], PALETTE[4], PALETTE[0]]
for k, (lab, vals) in enumerate(schemes.items()):
    bars = ax.bar(x + (k - 1.5) * w, vals, w * 0.92, label=lab, color=cols[k],
                  alpha=0.55 if k < 3 else 1.0)
    for b, val in zip(bars, vals):
        ax.text(b.get_x() + b.get_width() / 2, val + 1.2, f"{val:.0f}", ha="center", va="bottom",
                fontsize=8, color=INK)
ax.set_xticks(x)
ax.set_xticklabels(["firm A", "firm B", "firm C\n(worthless alone)"])
ax.set_ylabel("share of the joint value 90")
ax.set_ylim(0, 72)
ax.set_title("Who deserves what? Four ways to split v(ABC) = 90")
ax.legend(loc="upper left", fontsize=8.5, ncol=2)
ax.grid(axis="x", visible=False)
save(fig, "ch44-shapley-bars")

# ---------------------------------------------------------------- ranked-choice flow
fig, ax = plt.subplots(figsize=(8.6, 4.4), constrained_layout=True)
ax.axis("off")
ax.set_xlim(-0.3, 10.3)
ax.set_ylim(-8, 118)
C = {"A": PALETTE[0], "B": PALETTE[1], "C": PALETTE[2]}
x0, x1, bw = 1.0, 7.0, 1.4


def stack(xl, counts):
    """Draw a vertical stacked bar; return y-ranges per candidate."""
    y, out = 0, {}
    for c, n in counts:
        ax.add_patch(plt.Rectangle((xl, y), bw, n, fc=C[c], ec="white", lw=2))
        out[c] = (y, y + n)
        y += n + 2
    return out


r1 = stack(x0, [("A", 40), ("B", 35), ("C", 25)])
r2 = stack(x1, [("A", 45), ("B", 55)])
for c, (a, b) in r1.items():
    ax.text(x0 - 0.15, (a + b) / 2, f"{c}\n{b - a}", ha="right", va="center", fontsize=10, color=C[c],
            fontweight="bold")
for c, (a, b) in r2.items():
    ax.text(x1 + bw + 0.15, (a + b) / 2, f"{c}\n{b - a}", ha="left", va="center", fontsize=10, color=C[c],
            fontweight="bold")


def ribbon(ya, yb, color, alpha=0.25):
    t = np.linspace(0, 1, 60)
    s = 3 * t ** 2 - 2 * t ** 3
    xs = x0 + bw + t * (x1 - x0 - bw)
    lo = ya[0] + s * (yb[0] - ya[0])
    hi = ya[1] + s * (yb[1] - ya[1])
    ax.fill_between(xs, lo, hi, color=color, alpha=alpha, lw=0)


ribbon((0, 40), (0, 40), C["A"])                       # A keeps its 40
ribbon((42, 77), (47, 82), C["B"])                     # B keeps its 35
ribbon((79, 99), (82, 102), C["C"], alpha=0.45)        # 20 of C's ballots go to B
ribbon((99, 104), (40, 45), C["C"], alpha=0.45)        # 5 of C's ballots go to A
ax.text(x0 + bw / 2, -5, "Round 1\n(first choices)", ha="center", va="top", fontsize=9.5, color=INK)
ax.text(x1 + bw / 2, -5, "Round 2\n(C eliminated)", ha="center", va="top", fontsize=9.5, color=INK)
ax.text(4.75, 107, "C has the fewest first choices and is eliminated:\nits 25 ballots transfer, 20 to B and 5 to A", ha="center", va="bottom", fontsize=9,
        color=MUTED)
ax.text(x1 + bw / 2, 108, "B wins 55 : 45", ha="center", va="bottom", fontsize=10.5, fontweight="bold", color=C["B"])
ax.text(4.75, 18, "yet C beats A 60 : 40\nand beats B 65 : 35\nhead to head", ha="center", va="center",
        fontsize=9, color=C["C"])
save(fig, "ch44-rcv-flow")
