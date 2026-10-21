import sys, pathlib; sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
from style import plt, np, PALETTE, INK, MUTED, new_figure, diagram, node, box, arrow, save
from matplotlib.patches import Rectangle, Ellipse, Polygon, FancyArrowPatch


def self_loop(ax, x, y, r, text, angle=90, color=MUTED):
    """A self-loop drawn above (angle=90) or below (angle=-90) a node of radius r."""
    a = np.deg2rad(angle)
    p0 = (x + r * np.cos(a + 0.5), y + r * np.sin(a + 0.5))
    p1 = (x + r * np.cos(a - 0.5), y + r * np.sin(a - 0.5))
    ax.add_patch(FancyArrowPatch(p0, p1, connectionstyle="arc3,rad=-2.4", arrowstyle="-|>",
                                 mutation_scale=12, color=color, lw=1.6, shrinkA=0, shrinkB=0, zorder=2))
    ax.text(x + 1.25 * r * np.cos(a) * 1.9, y + r * np.sin(a) * 2.65, text, ha="center",
            va="center", fontsize=9, color=INK)


def edge(ax, p, q, text, curve=0.0, off=(0, 0.18), color=MUTED, r=0.62):
    d = np.subtract(q, p) / np.hypot(*np.subtract(q, p))
    t = -np.arctan(2 * curve)
    rot = np.array([[np.cos(t), -np.sin(t)], [np.sin(t), np.cos(t)]])
    p0 = np.add(p, r * (rot @ d)); q0 = np.subtract(q, r * (rot.T @ d))
    ax.add_patch(FancyArrowPatch(p0, q0, connectionstyle=f"arc3,rad={curve}", arrowstyle="-|>",
                                 mutation_scale=12, color=color, lw=1.6, shrinkA=0, shrinkB=0, zorder=2))
    mx, my = (p[0] + q[0]) / 2, (p[1] + q[1]) / 2
    ax.text(mx + off[0], my + off[1], text, ha="center", va="center", fontsize=9, color=INK)


# ---- Turing machine: tape, head, and the state diagram of a binary incrementer
fig, ax = diagram(figsize=(9.0, 5.6))
ax.set_xlim(0, 12); ax.set_ylim(0, 7.4)
tape = ["_", "_", "1", "0", "1", "1", "_", "_", "_"]
for k, s in enumerate(tape):
    ax.add_patch(Rectangle((1.5 + k, 5.6), 1, 0.9, fc="#f9fafb" if s == "_" else "#e5edff", ec=INK, lw=1.4))
    ax.text(2.0 + k, 6.05, s if s != "_" else "\u2423", ha="center", va="center", fontsize=13, color=INK)
ax.text(1.3, 6.05, "\u2026", ha="right", va="center", fontsize=14, color=MUTED)
ax.text(10.7, 6.05, "\u2026", ha="left", va="center", fontsize=14, color=MUTED)
hx = 2.0 + 5
ax.add_patch(Polygon([[hx, 5.55], [hx - 0.35, 4.95], [hx + 0.35, 4.95]], fc=PALETTE[1], ec=PALETTE[1]))
ax.text(hx, 4.7, "head, in state carry", ha="center", va="top", fontsize=9, color=PALETTE[1])
ax.text(6.0, 7.05, "infinite tape (input 1011 = 11; it will become 1100 = 12)", ha="center", fontsize=9.5, color=MUTED)
S = {"right": (2.2, 2.0), "carry": (6.0, 2.0), "halt": (9.8, 2.0)}
for name, (x, y) in S.items():
    node(ax, x, y, name, color=PALETTE[3] if name == "halt" else PALETTE[0], r=0.62, fontsize=10,
         lw=3 if name == "halt" else 2)
ax.annotate("", xy=(S["right"][0] - 0.62, 2.0), xytext=(0.5, 2.0),
            arrowprops=dict(arrowstyle="-|>", color=INK, lw=1.6))
ax.text(0.55, 2.25, "start", fontsize=9, color=INK)
self_loop(ax, *S["right"], 0.62, "0 \u2192 0, R\n1 \u2192 1, R", angle=90)
self_loop(ax, *S["carry"], 0.62, "1 \u2192 0, L", angle=90)
edge(ax, S["right"], S["carry"], "\u2423 \u2192 \u2423, L", off=(0, 0.25))
edge(ax, S["carry"], S["halt"], "0 \u2192 1, L\n\u2423 \u2192 1, L", off=(0, 0.42))
ax.text(6.0, 0.45, "each arrow: read symbol \u2192 write symbol, move head (L/R)", ha="center", fontsize=9, color=MUTED)
save(fig, "ch30-turing-machine")

# ---- DFA: binary strings whose value is divisible by 3
fig, ax = diagram(figsize=(8.0, 4.0))
ax.set_xlim(0, 10); ax.set_ylim(0, 5)
Q = {0: (2.0, 2.2), 1: (5.0, 2.2), 2: (8.0, 2.2)}
for r, (x, y) in Q.items():
    node(ax, x, y, f"$r={r}$", color=PALETTE[2] if r == 0 else PALETTE[0], r=0.55, fontsize=10.5)
ax.add_patch(plt.Circle(Q[0], 0.68, fc="none", ec=PALETTE[2], lw=1.8, zorder=3))
ax.annotate("", xy=(Q[0][0] - 0.68, 2.2), xytext=(0.3, 2.2), arrowprops=dict(arrowstyle="-|>", color=INK, lw=1.6))
ax.text(0.3, 2.45, "start", fontsize=9, color=INK)
self_loop(ax, *Q[0], 0.68, "0", angle=90)
self_loop(ax, *Q[2], 0.55, "1", angle=90)
edge(ax, Q[0], Q[1], "1", curve=-0.35, off=(0, 0.75), r=0.6)
edge(ax, Q[1], Q[0], "1", curve=-0.35, off=(0, -0.75), r=0.6)
edge(ax, Q[1], Q[2], "0", curve=-0.35, off=(0, 0.75), r=0.55)
edge(ax, Q[2], Q[1], "0", curve=-0.35, off=(0, -0.75), r=0.55)
ax.text(5.0, 0.35, r"state = (value so far) mod 3;  reading bit $b$: $r \to (2r+b)\ \mathrm{mod}\ 3$;  accept in $r=0$",
        ha="center", fontsize=9, color=MUTED)
save(fig, "ch30-dfa")

# ---- complexity map
fig, ax = diagram(figsize=(9.0, 5.6))
ax.set_xlim(0, 12); ax.set_ylim(0, 7.4)
ax.add_patch(Rectangle((0.2, 0.2), 11.6, 7.0, fc="#fef2f2", ec=PALETTE[1], lw=1.5))
ax.text(0.5, 6.85, "all decision problems (uncountably many)", fontsize=10, color=PALETTE[1], weight="bold")
ax.text(10.0, 0.55, "undecidable:\nHALT, Rice-type\nquestions", ha="center", fontsize=9, color=PALETTE[1])
ax.add_patch(Ellipse((5.0, 3.6), 9.0, 5.8, fc="#f5f3ff", ec=PALETTE[3], lw=1.8))
ax.text(5.0, 6.15, "decidable", ha="center", fontsize=10, color=PALETTE[3], weight="bold")
ax.add_patch(Ellipse((4.6, 3.3), 7.2, 4.4, fc="#ecfeff", ec=PALETTE[5], lw=1.8))
ax.text(4.6, 5.15, "EXP", ha="center", fontsize=10, color=PALETTE[5], weight="bold")
ax.add_patch(Ellipse((4.2, 3.0), 5.2, 3.2, fc="#fff7ed", ec=PALETTE[4], lw=1.8))
ax.text(5.4, 4.15, "NP", ha="center", fontsize=10, color=PALETTE[4], weight="bold")
ax.add_patch(Ellipse((3.3, 2.75), 2.4, 1.8, fc="#eff6ff", ec=PALETTE[0], lw=1.8))
ax.text(3.3, 3.15, "P", ha="center", fontsize=11, color=PALETTE[0], weight="bold")
ax.text(3.3, 2.45, "sorting\nshortest path\nprimality", ha="center", fontsize=7.5, color=INK)
ax.add_patch(Ellipse((5.75, 2.75), 1.7, 1.5, fc="#ffedd5", ec=PALETTE[4], lw=1.5, ls="--"))
ax.text(5.75, 2.75, "NP-complete\nSAT, 3-colour,\nclique, TSP", ha="center", va="center", fontsize=7.5, color=INK)
ax.text(7.5, 3.0, "n×n chess,\nGo:\nEXP-\ncomplete", ha="center", va="center", fontsize=7.5, color=INK)
ax.text(4.2, 1.6, "factoring?  graph isomorphism?", ha="center", fontsize=7.5, color=MUTED)
ax.text(11.6, 6.5, "picture assumes P \u2260 NP;\nif P = NP the two\ninner ovals merge", ha="right", va="top", fontsize=8.5, color=MUTED)
save(fig, "ch30-complexity-map")
