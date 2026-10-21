import sys, pathlib; sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
from style import plt, np, PALETTE, INK, MUTED, new_figure, diagram, node, box, arrow, save
from matplotlib.patches import FancyArrowPatch

# ---- Hoare-logic flow through the division loop, invariant highlighted
fig, ax = diagram(figsize=(9.6, 6.6))
ax.set_xlim(0, 13); ax.set_ylim(0, 9.0)
BLUE, RED, GREEN, PURPLE = PALETTE[0], PALETTE[1], PALETTE[2], PALETTE[3]
box(ax, 4.0, 8.3, r"$\{P\}:\ x\geq 0\ \wedge\ y>0$", w=4.6, h=0.7, color=MUTED, fc="#f9fafb", fontsize=10.5)
box(ax, 4.0, 7.0, "q := 0;   r := x", w=3.6, h=0.7, color=INK, fontsize=10.5)
box(ax, 4.0, 5.4, r"invariant $I$:  $x = q\,y + r\ \wedge\ r\geq 0$", w=6.4, h=0.95, color=PALETTE[4],
    fc="#fff7ed", fontsize=11.5, lw=3)
box(ax, 4.0, 3.8, r"guard:  $r\geq y$ ?", w=3.4, h=0.7, color=INK, fontsize=10.5)
box(ax, 10.0, 3.8, "r := r - y;   q := q + 1", w=4.2, h=0.7, color=INK, fontsize=10.5)
box(ax, 4.0, 1.6, r"$\{Q\}:\ x = q\,y + r\ \wedge\ 0\leq r<y$", w=5.4, h=0.7, color=MUTED, fc="#f9fafb",
    fontsize=10.5)
arrow(ax, 4.0, 8.3, 4.0, 7.0, color=INK, shrink=16)
arrow(ax, 4.0, 7.0, 4.0, 5.4, color=INK, shrink=18)
arrow(ax, 4.0, 5.4, 4.0, 3.8, color=INK, shrink=20)
arrow(ax, 4.0, 3.8, 4.0, 1.6, color=INK, shrink=18)
ax.text(4.15, 2.75, "no", fontsize=9.5, color=INK)
ax.annotate("", xy=(7.9, 3.8), xytext=(5.7, 3.8), arrowprops=dict(arrowstyle="-|>", color=INK, lw=1.6))
ax.text(6.8, 3.95, "yes", ha="center", fontsize=9.5, color=INK)
ax.add_patch(FancyArrowPatch((10.0, 4.15), (7.2, 5.5), connectionstyle="arc3,rad=0.35", arrowstyle="-|>",
                             mutation_scale=14, color=INK, lw=1.6))
# verification conditions
vc = dict(ha="left", va="center", fontsize=9.2)
ax.text(6.4, 7.0, r"VC1:  $P \Rightarrow x = 0\cdot y + x\ \wedge\ x\geq 0$", color=BLUE, **vc)
ax.text(8.3, 5.95, "VC2:  $I \\wedge r\\geq y \\Rightarrow$\n$x=(q+1)\\,y+(r-y)\\ \\wedge\\ r-y\\geq 0$", color=BLUE, **vc)
ax.text(7.1, 1.6, r"VC3:  $I \wedge r<y \Rightarrow Q$", color=BLUE, **vc)
ax.text(10.0, 2.8, r"variant $r$ drops by $y>0$" + "\neach pass, and stays $\\geq 0$", ha="center", va="center",
        fontsize=9, color=PURPLE)
ax.text(0.2, 0.45, "Three implications of ordinary arithmetic replace infinitely many test runs.", fontsize=9.5,
        color=MUTED)
save(fig, "ch67-hoare-loop")

# ---- model-checking state graph of the naive lock, bad state reachable
fig, ax = diagram(figsize=(10.0, 5.2))
ax.set_xlim(0, 13); ax.set_ylim(-0.5, 6.6)
name = {0: "I", 1: "T", 2: "C"}
depth_states = [[(0, 0, 0)], [(1, 0, 0), (0, 1, 0)], [(2, 0, 1), (1, 1, 0), (0, 2, 1)],
                [(2, 1, 1), (1, 2, 1)], [(2, 2, 1)]]
pos = {}
for d, layer in enumerate(depth_states):
    ys = np.linspace(3.3 + 1.25 * (len(layer) - 1), 3.3 - 1.25 * (len(layer) - 1), len(layer))
    for s, y in zip(layer, ys):
        pos[s] = (1.2 + 2.6 * d, y)
cex = [(0, 0, 0), (1, 0, 0), (1, 1, 0), (2, 1, 1), (2, 2, 1)]
fwd = [((0, 0, 0), (1, 0, 0)), ((0, 0, 0), (0, 1, 0)), ((1, 0, 0), (2, 0, 1)), ((1, 0, 0), (1, 1, 0)),
       ((0, 1, 0), (1, 1, 0)), ((0, 1, 0), (0, 2, 1)), ((1, 1, 0), (2, 1, 1)), ((1, 1, 0), (1, 2, 1)),
       ((2, 1, 1), (2, 2, 1)), ((1, 2, 1), (2, 2, 1))]
back = [((2, 0, 1), (0, 0, 0)), ((0, 2, 1), (0, 0, 0)), ((2, 1, 1), (0, 1, 0)), ((1, 2, 1), (1, 0, 0))]
R = 0.58
for s, t in fwd:
    on = s in cex and t in cex and cex.index(t) == cex.index(s) + 1
    ax.add_patch(FancyArrowPatch(pos[s], pos[t], arrowstyle="-|>", mutation_scale=12, shrinkA=R * 72 / 1.3,
                                 shrinkB=R * 72 / 1.3, color=PALETTE[1] if on else MUTED, lw=2.6 if on else 1.3,
                                 zorder=2))
for s, t in back:
    ax.add_patch(FancyArrowPatch(pos[s], pos[t], arrowstyle="-|>", mutation_scale=10, shrinkA=R * 72 / 1.3,
                                 shrinkB=R * 72 / 1.3, color="#c4c8cf", lw=1.0, ls="--",
                                 connectionstyle="arc3,rad=0.45", zorder=1))
for s, (x, y) in pos.items():
    bad = s == (2, 2, 1)
    node(ax, x, y, f"{name[s[0]]} {name[s[1]]}\nlock={s[2]}", color=PALETTE[1] if bad else PALETTE[0],
         r=R, fc="#fee2e2" if bad else "white", fontsize=8.5, lw=3 if bad else 2)
ax.text(pos[(2, 2, 1)][0], pos[(2, 2, 1)][1] - 0.95, "BAD: both in\ncritical section", ha="center", va="top",
        fontsize=9, color=PALETTE[1], weight="bold")
ax.annotate("", xy=(pos[(0, 0, 0)][0] - R, 3.3), xytext=(0.1, 3.3), arrowprops=dict(arrowstyle="-|>", color=INK))
ax.text(0.1, 3.5, "start", fontsize=8.5, color=INK)
ax.text(6.5, -0.35, "I = idle, T = has seen lock free, C = critical.   Red: the shortest counterexample found by "
        "breadth-first search.   Dashed: release steps.", ha="center", fontsize=8.3, color=MUTED)
save(fig, "ch67-state-graph")
