"""Chapter 39 figures: payoff matrices with best responses marked; best-response curves of the penalty kick."""
import sys, pathlib; sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
from style import plt, np, PALETTE, INK, MUTED, new_figure, diagram, node, box, arrow, save
from matplotlib.patches import Rectangle, Ellipse

BLUE, RED, GREEN = PALETTE[0], PALETTE[1], PALETTE[2]


def matrix(ax, x0, y0, rows, cols, A, B, title, cell=1.3, ne=()):
    """Draw a bimatrix with its top-left corner at (x0, y0); circle best responses."""
    m, n = A.shape
    ax.text(x0 + n * cell / 2, y0 + 0.95, title, ha="center", va="bottom", fontsize=11, weight="bold", color=INK)
    for j, c in enumerate(cols):
        ax.text(x0 + (j + 0.5) * cell, y0 + 0.3, c, ha="center", va="center", fontsize=10, color=RED, weight="bold")
    for i, r in enumerate(rows):
        ax.text(x0 - 0.35, y0 - (i + 0.5) * cell, r, ha="center", va="center", fontsize=10, color=BLUE, weight="bold")
    for i in range(m):
        for j in range(n):
            x, y = x0 + j * cell, y0 - (i + 1) * cell
            fc = "#dcfce7" if (i, j) in ne else "white"
            ax.add_patch(Rectangle((x, y), cell, cell, fc=fc, ec=INK, lw=1.2))
            xa, xb, yc = x + 0.28 * cell, x + 0.72 * cell, y + 0.5 * cell
            ax.text(xa, yc, f"{A[i, j]:g}", ha="center", va="center", fontsize=10.5, color=BLUE)
            ax.text(x + 0.5 * cell, yc, ",", ha="center", va="center", fontsize=10.5, color=INK)
            ax.text(xb, yc, f"{B[i, j]:g}", ha="center", va="center", fontsize=10.5, color=RED)
            if A[i, j] == A[:, j].max():
                ax.add_patch(Ellipse((xa, yc), 0.36 * cell, 0.36 * cell, fc="none", ec=BLUE, lw=1.4))
            if B[i, j] == B[i, :].max():
                ax.add_patch(Ellipse((xb, yc), 0.36 * cell, 0.36 * cell, fc="none", ec=RED, lw=1.4))


fig, ax = plt.subplots(figsize=(9.6, 4.4), constrained_layout=True)
ax.set_xlim(-0.6, 10.6); ax.set_ylim(-4.6, 1.6); ax.set_aspect("equal"); ax.axis("off")
A = np.array([[4, 1, 3], [2, 2, 1], [3, 0, 5]]); B = np.array([[3, 1, 0], [1, 4, 2], [0, 2, 1]])
matrix(ax, 0.4, 0, ["T", "M", "B"], ["L", "C", "R"], A, B, "(a) Two pure equilibria", ne={(0, 0), (1, 1)})
S = np.array([[58, 95], [93, 70]])
matrix(ax, 6.6, 0, ["kick L", "kick R"], ["dive L", "dive R"], S, 100 - S, "(b) Penalty kick: no pure equilibrium", cell=1.6)
ax.text(6.6 - 1.05, -0.8, "", fontsize=9)
# cycle arrows in (b)
c = 1.6; x0 = 6.6
pts = [(x0 + 0.5 * c, -0.5 * c), (x0 + 1.5 * c, -0.5 * c), (x0 + 1.5 * c, -1.5 * c), (x0 + 0.5 * c, -1.5 * c)]
# keeper moves along rows, kicker along columns: from (L,L) kicker -> R ; etc.
seq = [(0, 3), (3, 2), (2, 1), (1, 0)]
for a, b in seq:
    (xa, ya), (xb, yb) = pts[a], pts[b]
    ax.annotate("", xy=(xb, yb), xytext=(xa, ya), arrowprops=dict(arrowstyle="-|>", color=GREEN, lw=1.6,
                shrinkA=22, shrinkB=22, connectionstyle="arc3,rad=0.25"), zorder=5)
ax.text(0.4 + 1.95, -4.35, "circled = best response; shaded = Nash equilibrium", ha="center", fontsize=9, color=MUTED)
ax.text(x0 + c, -3.75, "scoring chance (%), kicker , keeper\narrows: someone always wants to switch", ha="center",
        va="top", fontsize=9, color=MUTED)
save(fig, "ch39-payoff-matrices")

# best-response curves
fig, ax = new_figure(figsize=(5.6, 5.0))
ps, qs = 23 / 60, 5 / 12
ax.plot([1, 1, 0, 0], [0, qs, qs, 1], color=BLUE, lw=2.6, label="kicker's best response $p^*(q)$")
ax.plot([0, ps, ps, 1], [0, 0, 1, 1], color=RED, lw=2.6, ls="--", label="keeper's best response $q^*(p)$")
ax.plot(ps, qs, "o", ms=10, color=GREEN, zorder=5)
ax.annotate(f"Nash equilibrium\n$p={ps:.3f}$, $q={qs:.3f}$", xy=(ps, qs), xytext=(0.5, 0.18), fontsize=10,
            arrowprops=dict(arrowstyle="-|>", color=INK, lw=1.2))
ax.set_xlim(-0.03, 1.03); ax.set_ylim(-0.03, 1.03)
ax.set_xlabel("$p$ = probability the kicker shoots left")
ax.set_ylabel("$q$ = probability the keeper dives left")
ax.legend(loc="upper right", bbox_to_anchor=(1.0, 0.94), fontsize=9)
ax.set_title("The only point where both are best responding")
save(fig, "ch39-best-response")
