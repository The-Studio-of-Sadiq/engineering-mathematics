import sys, pathlib; sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
from style import plt, np, PALETTE, INK, MUTED, new_figure, diagram, node, box, arrow, save
from matplotlib.patches import Rectangle

ON, OFF = "#bbf7d0", "#f3f4f6"


def grid(ax, table, title, rows, cols, rlab, clab):
    """Draw a coloured truth grid: table[i][j] in {0,1}."""
    nr, nc = len(table), len(table[0])
    for i in range(nr):
        for j in range(nc):
            v = table[i][j]
            ax.add_patch(Rectangle((j, nr - 1 - i), 1, 1, fc=ON if v else OFF, ec=INK, lw=1.5))
            ax.text(j + 0.5, nr - 1 - i + 0.5, str(v), ha="center", va="center", fontsize=18,
                    weight="bold", color=PALETTE[2] if v else MUTED)
    for j, c in enumerate(cols):
        ax.text(j + 0.5, nr + 0.18, c, ha="center", va="bottom", fontsize=10, color=INK)
    for i, r in enumerate(rows):
        ax.text(-0.15, nr - 1 - i + 0.5, r, ha="right", va="center", fontsize=10, color=INK)
    if clab:
        ax.text(nc / 2, nr + 0.62, clab, ha="center", va="bottom", fontsize=10, color=MUTED)
    if rlab:
        ax.text(-0.75, nr / 2, rlab, ha="center", va="center", fontsize=10, color=MUTED, rotation=90)
    ax.set_title(title, pad=6)
    ax.set_xlim(-1.0, 2.1)
    ax.set_ylim(-0.1, 2.75)
    ax.set_aspect("equal")
    ax.axis("off")


fig, axs = new_figure(1, 4, figsize=(10.0, 3.2))
grid(axs[0], [[1], [0]], r"NOT $\neg x$", ["x=0", "x=1"], [r"$\neg x$"], None, None)
grid(axs[1], [[0, 0], [0, 1]], r"AND $x\wedge y$", ["x=0", "x=1"], ["y=0", "y=1"], None, None)
grid(axs[2], [[0, 1], [1, 1]], r"OR $x\vee y$", ["x=0", "x=1"], ["y=0", "y=1"], None, None)
grid(axs[3], [[0, 1], [1, 0]], r"XOR $x\oplus y$", ["x=0", "x=1"], ["y=0", "y=1"], None, None)
save(fig, "ch29-truth-tables")

# ---- loop-invariant diagram for binary search
fig, ax = diagram(figsize=(9.0, 5.0))
ax.set_xlim(0, 12); ax.set_ylim(0, 6.6)
vals = [2, 5, 8, 12, 16, 23, 38, 56, 72, 91]
lo, hi = 3, 7
x0, w = 1.0, 1.0
for k, v in enumerate(vals):
    if k < lo:
        fc = "#fee2e2"
    elif k >= hi:
        fc = "#dbeafe"
    else:
        fc = "#fef9c3"
    ax.add_patch(Rectangle((x0 + k * w, 4.3), w, 0.9, fc=fc, ec=INK, lw=1.4))
    ax.text(x0 + k * w + 0.5, 4.75, str(v), ha="center", va="center", fontsize=11, color=INK)
    ax.text(x0 + k * w + 0.5, 4.05, str(k), ha="center", va="center", fontsize=8, color=MUTED)
ax.text(0.85, 4.05, "index", ha="right", va="center", fontsize=8, color=MUTED)
ax.annotate("lo", xy=(x0 + lo * w, 5.25), xytext=(x0 + lo * w, 5.85), ha="center", fontsize=10,
            color=INK, arrowprops=dict(arrowstyle="-|>", color=INK))
ax.annotate("hi", xy=(x0 + hi * w, 5.25), xytext=(x0 + hi * w, 5.85), ha="center", fontsize=10,
            color=INK, arrowprops=dict(arrowstyle="-|>", color=INK))
ax.text(x0 + lo * w / 2, 5.45, "all < target", ha="center", fontsize=9, color=PALETTE[1])
ax.text(x0 + (hi + 10) * w / 2, 5.45, "all > target", ha="center", fontsize=9, color=PALETTE[0])
ax.text(x0 + (lo + hi) * w / 2, 3.65, "target, if present, is in here", ha="center", fontsize=9, color="#a16207")
box(ax, 2.0, 1.7, "Initialise\nlo = 0, hi = n\ninvariant holds:\nnothing excluded", w=3.0, h=1.7, color=PALETTE[2])
box(ax, 6.0, 1.7, "Maintain\none step keeps it true;\nvariant hi - lo\nstrictly shrinks", w=3.0, h=1.7, color=PALETTE[4])
box(ax, 10.0, 1.7, "Terminate\nlo = hi: empty window\n+ invariant\n= correct answer", w=3.0, h=1.7, color=PALETTE[3])
arrow(ax, 3.5, 1.7, 4.5, 1.7, color=INK, shrink=2)
arrow(ax, 7.5, 1.7, 8.5, 1.7, color=INK, shrink=2)
ax.annotate("", xy=(6.6, 2.58), xytext=(5.4, 2.58), arrowprops=dict(arrowstyle="-|>", color=PALETTE[4], lw=1.6,
            connectionstyle="arc3,rad=-1.2"))
ax.text(7.05, 2.95, "repeat", ha="center", fontsize=8.5, color=PALETTE[4])
save(fig, "ch29-loop-invariant")

# ---- Hasse diagrams: divisors of 30 and subsets of {a,b,c}
fig, axs = new_figure(1, 2, figsize=(9.0, 4.2))
pos = {0: (0, 0), 1: (-1.4, 1.2), 2: (0, 1.2), 4: (1.4, 1.2), 3: (-1.4, 2.4), 5: (0, 2.4), 6: (1.4, 2.4), 7: (0, 3.6)}
labels_div = {0: "1", 1: "2", 2: "3", 4: "5", 3: "6", 5: "10", 6: "15", 7: "30"}
names = "abc"
labels_set = {m: "{" + ",".join(names[i] for i in range(3) if m >> i & 1) + "}" if m else r"$\varnothing$" for m in range(8)}
for ax, labels, title, col in [(axs[0], labels_div, "divisors of 30, ordered by 'divides'", PALETTE[0]),
                               (axs[1], labels_set, r"subsets of {a,b,c}, ordered by $\subseteq$", PALETTE[3])]:
    for m in range(8):
        for i in range(3):
            if not m >> i & 1:
                p, q = pos[m], pos[m | 1 << i]
                ax.plot([p[0], q[0]], [p[1], q[1]], color=MUTED, lw=1.5, zorder=1)
    for m, (x, y) in pos.items():
        ax.add_patch(plt.Circle((x, y), 0.38, fc="white", ec=col, lw=2, zorder=3))
        ax.text(x, y, labels[m], ha="center", va="center", fontsize=9.5 if ax is axs[0] else 8, color=INK, zorder=4)
    ax.set_xlim(-2.1, 2.1); ax.set_ylim(-0.6, 4.2); ax.set_aspect("equal"); ax.axis("off")
    ax.set_title(title, fontsize=10.5)
axs[0].text(0, -0.55, "meet = gcd, join = lcm", ha="center", fontsize=9, color=MUTED)
axs[1].text(0, -0.55, r"meet = $\cap$, join = $\cup$", ha="center", fontsize=9, color=MUTED)
save(fig, "ch29-hasse")
