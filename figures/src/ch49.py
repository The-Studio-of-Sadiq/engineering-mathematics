"""Figures for Chapter 49 — Integer Programming and Combinatorial Optimization."""
import sys, pathlib; sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
from style import plt, np, PALETTE, INK, MUTED, new_figure, diagram, node, box, arrow, save

BLUE, RED, GREEN, PURPLE, ORANGE = PALETTE[:5]

# ------------------------------------------- figure 1: LP relaxation polytope and lattice points
# max 5x + 8y  s.t.  x + y <= 6,  5x + 9y <= 45,  x, y >= 0 integer
fig, ax = new_figure(figsize=(7.2, 5.2))
lp = np.array([[0, 0], [6, 0], [2.25, 3.75], [0, 5]])
ax.fill(lp[:, 0], lp[:, 1], color=BLUE, alpha=0.12, label="LP relaxation")
ax.plot(*np.r_[lp, lp[:1]].T, color=BLUE, lw=2)
hull = np.array([[0, 0], [6, 0], [3, 3], [0, 5]])
ax.plot(*np.r_[hull, hull[:1]].T, color=GREEN, lw=2, ls="--", label="convex hull of integer points")
for x in range(8):
    for y in range(7):
        ok = x + y <= 6 and 5 * x + 9 * y <= 45
        ax.plot(x, y, "o", ms=6 if ok else 4, color=INK if ok else "#cbd5e1", zorder=4)
xs = np.linspace(-0.5, 7.5, 2)
ax.plot(xs, (15.75 - 2 * xs) / 3, color=PURPLE, lw=1.0, ls=":")
ax.plot(xs, (15 - 2 * xs) / 3, color=PURPLE, lw=1.6, label=r"cut $2x+3y\leq15$")
ax.plot(xs, (41.25 - 5 * xs) / 8, color=RED, lw=1.1, ls="--")
ax.plot(2.25, 3.75, "*", ms=15, color=RED, zorder=6)
ax.annotate("LP optimum (2.25, 3.75)\nvalue 41.25", xy=(2.25, 3.75), xytext=(3.6, 5.55), fontsize=9.5,
            arrowprops=dict(arrowstyle="-", color=INK, lw=0.8))
ax.plot(0, 5, "s", ms=10, mfc="none", mec=GREEN, mew=2.2, zorder=6)
ax.annotate("integer optimum (0, 5)\nvalue 40", xy=(0, 5), xytext=(0.3, 6.15), fontsize=9.5,
            arrowprops=dict(arrowstyle="-", color=INK, lw=0.8))
ax.plot(2, 4, "x", ms=10, mew=2.2, color=ORANGE, zorder=6)
ax.annotate("rounded (2, 4):\ninfeasible", xy=(2, 4), xytext=(3.1, 4.45), fontsize=9, color=ORANGE,
            arrowprops=dict(arrowstyle="-", color=ORANGE, lw=0.8))
ax.plot(2, 3, "x", ms=10, mew=2.2, color=ORANGE, zorder=6)
ax.annotate("rounded down (2, 3):\nvalue 34", xy=(2, 3), xytext=(4.4, 3.3), fontsize=9, color=ORANGE,
            arrowprops=dict(arrowstyle="-", color=ORANGE, lw=0.8))
ax.set_xlim(-0.5, 7.5); ax.set_ylim(-0.5, 6.9); ax.set_aspect("equal")
ax.set_xlabel("$x$"); ax.set_ylabel("$y$")
ax.set_title(r"maximise $5x+8y$:  $x+y\leq6$,  $5x+9y\leq45$,  $x,y\in\mathbb{Z}_{\geq0}$", fontsize=11)
ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.11), ncol=3, fontsize=8.5)
save(fig, "ch49-relaxation")

# ------------------------------------------- figure 2: the branch-and-bound tree for the same problem
fig, ax = diagram(figsize=(9.5, 6.0))
ax.set_xlim(0, 12); ax.set_ylim(-0.2, 7.5)
W, H = 2.35, 0.95
nodes = {
    "R":    (6.0, 6.9, "root\n(2.25, 3.75)  41.25", BLUE),
    "L":    (3.6, 5.2, "(2, 3.89)  41.11", BLUE),
    "Rt":   (9.3, 5.2, "(3, 3)  39\nsolved last: 39 < 40, pruned", MUTED),
    "LL":   (1.5, 3.5, "(2, 3)  34\ninteger: incumbent", ORANGE),
    "LR":   (5.4, 3.5, "(1.8, 4)  41", BLUE),
    "LRL":  (3.6, 1.8, "(1, 4.44)  40.56", BLUE),
    "LRR":  (7.6, 1.8, "infeasible", MUTED),
    "LRLL": (1.5, 0.45, "(1, 4)  37\nworse: pruned", MUTED),
    "LRLR": (5.6, 0.45, "(0, 5)  40\nnew incumbent", GREEN),
}
edges = [("R", "L", r"$x\leq2$"), ("R", "Rt", r"$x\geq3$"), ("L", "LL", r"$y\leq3$"), ("L", "LR", r"$y\geq4$"),
         ("LR", "LRL", r"$x\leq1$"), ("LR", "LRR", r"$x\geq2$"), ("LRL", "LRLL", r"$y\leq4$"),
         ("LRL", "LRLR", r"$y\geq5$")]
for a, b, lab in edges:
    xa, ya = nodes[a][:2]; xb, yb = nodes[b][:2]
    ax.plot([xa, xb], [ya - H / 2, yb + H / 2], color=MUTED, lw=1.3, zorder=1)
    ax.text((xa + xb) / 2 + (0.12 if xb > xa else -0.12), (ya + yb) / 2, lab, fontsize=9,
            ha="left" if xb > xa else "right", va="center", color=INK)
for key, (x, y, txt, col) in nodes.items():
    box(ax, x, y, txt, w=W, h=H, color=col, fontsize=8.5)
save(fig, "ch49-branch-and-bound")
