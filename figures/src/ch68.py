import sys, pathlib; sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
from style import plt, np, PALETTE, INK, MUTED, new_figure, diagram, node, box, arrow, save
from matplotlib.patches import Rectangle, Ellipse, FancyArrowPatch
import math

# ---- Russell's set as a diagonal: membership table, flipped diagonal, self-reference loop
fig, ax = diagram(figsize=(9.6, 5.2))
ax.set_xlim(0, 13); ax.set_ylim(0, 6.8)
M = [[0, 1, 0, 0, 1], [1, 1, 0, 1, 0], [0, 0, 0, 1, 1], [1, 0, 1, 1, 0], [0, 1, 1, 0, 0]]
n, x0, y0, c = 5, 1.6, 1.2, 0.9
for i in range(n):
    for j in range(n):
        diag = i == j
        ax.add_patch(Rectangle((x0 + j * c, y0 + (n - 1 - i) * c), c, c, fc="#fef3c7" if diag else "white",
                               ec=INK, lw=1.2))
        ax.text(x0 + j * c + c / 2, y0 + (n - 1 - i) * c + c / 2, "∈" if M[i][j] else "∉",
                ha="center", va="center", fontsize=13, color=PALETTE[4] if diag else INK)
    ax.text(x0 - 0.15, y0 + (n - 1 - i) * c + c / 2, f"$x_{i + 1}$", ha="right", va="center", fontsize=11)
    ax.text(x0 + i * c + c / 2, y0 + n * c + 0.15, f"$x_{i + 1}$", ha="center", va="bottom", fontsize=11)
ax.text(x0 + n * c / 2, y0 + n * c + 0.75, r"is the row set a member of the column set?", ha="center", fontsize=9.5,
        color=MUTED)
ax.text(x0 - 0.9, y0 + n * c / 2, "every set", rotation=90, ha="center", va="center", fontsize=9.5, color=MUTED)
# R row: the flipped diagonal
ry = y0 - 1.0
ax.text(x0 - 0.15, ry + c / 2, "$R$", ha="right", va="center", fontsize=12, color=PALETTE[1], weight="bold")
for j in range(n):
    ax.add_patch(Rectangle((x0 + j * c, ry), c, c, fc="#fee2e2", ec=PALETTE[1], lw=1.4))
    ax.text(x0 + j * c + c / 2, ry + c / 2, "∉" if M[j][j] else "∈", ha="center", va="center",
            fontsize=13, color=PALETTE[1])
ax.text(x0 + n * c + 0.2, ry + c / 2, "flip the diagonal: $x_j\\in R\\ \\Leftrightarrow\\ x_j\\notin x_j$", ha="left",
        va="center", fontsize=9.5, color=PALETTE[1])
# right panel: the self-reference loop
cx, cy = 10.2, 4.4
box(ax, cx, cy + 1.5, r"suppose $R\in R$", w=3.0, h=0.7, color=PALETTE[0], fontsize=10.5)
box(ax, cx, cy - 1.5, r"suppose $R\notin R$", w=3.0, h=0.7, color=PALETTE[0], fontsize=10.5)
ax.add_patch(FancyArrowPatch((cx + 1.55, cy + 1.3), (cx + 1.55, cy - 1.3), connectionstyle="arc3,rad=-0.5",
                             arrowstyle="-|>", mutation_scale=14, color=PALETTE[1], lw=1.8))
ax.add_patch(FancyArrowPatch((cx - 1.55, cy - 1.3), (cx - 1.55, cy + 1.3), connectionstyle="arc3,rad=-0.5",
                             arrowstyle="-|>", mutation_scale=14, color=PALETTE[1], lw=1.8))
ax.text(cx + 2.35, cy, "then by\ndefinition\n$R\\notin R$", ha="center", va="center", fontsize=9, color=INK)
ax.text(cx - 2.35, cy, "then by\ndefinition\n$R\\in R$", ha="center", va="center", fontsize=9, color=INK)
ax.text(cx, cy, "if $R$ is one of\nthe rows, its own\ndiagonal cell\ncontradicts itself", ha="center", va="center",
        fontsize=8.8, color=MUTED)
save(fig, "ch68-russell-diagonal")

# ---- map of what Godel's theorems say about a formal system
fig, ax = diagram(figsize=(10.0, 5.4))
ax.set_xlim(0, 13); ax.set_ylim(0.6, 7.3)
ax.add_patch(Rectangle((0.3, 0.9), 8.4, 6.2, fc="#f9fafb", ec=INK, lw=1.4))
ax.text(0.5, 6.8, "all sentences of arithmetic", fontsize=10, color=INK, weight="bold", va="center")
ax.add_patch(Rectangle((0.3, 0.9), 4.2, 5.6, fc="#ecfdf5", ec=PALETTE[2], lw=1.6))
ax.text(2.4, 6.2, "TRUE in $\\mathbb{N}$", ha="center", fontsize=10.5, color=PALETTE[2], weight="bold")
ax.add_patch(Rectangle((4.5, 0.9), 4.2, 5.6, fc="#fef2f2", ec=PALETTE[1], lw=1.6))
ax.text(6.6, 6.2, "FALSE in $\\mathbb{N}$", ha="center", fontsize=10.5, color=PALETTE[1], weight="bold")
ax.add_patch(Ellipse((2.2, 3.9), 3.2, 3.0, fc="#dbeafe", ec=PALETTE[0], lw=1.8))
ax.text(2.2, 4.5, "provable in $T$", ha="center", fontsize=10, color=PALETTE[0], weight="bold")
ax.text(2.2, 3.6, "$2+2=4$\n$\\forall n\\ (n+0=n)$\ninfinitely many primes", ha="center", va="center", fontsize=8.5,
        color=INK)
ax.add_patch(Ellipse((6.8, 3.9), 3.2, 3.0, fc="#ffe4e6", ec=PALETTE[1], lw=1.8, ls="--"))
ax.text(6.8, 4.5, "refutable in $T$", ha="center", fontsize=10, color=PALETTE[1], weight="bold")
ax.text(6.8, 3.6, "$2+2=5$\n$\\exists n\\ (n<0)$", ha="center", va="center", fontsize=8.5, color=INK)
node(ax, 3.4, 2.0, "$G$", color=PALETTE[4], r=0.38, fc="#fff7ed", fontsize=11, lw=2.4)
node(ax, 1.3, 2.0, "Con$(T)$", color=PALETTE[3], r=0.55, fc="#f5f3ff", fontsize=8.5, lw=2.4)
ax.text(2.35, 1.0, "true, yet neither provable nor refutable", ha="center", va="bottom", fontsize=7.8,
        color=MUTED)
# legend panel on the right
tx = 9.0
ax.text(tx, 6.8, "If $T$ is consistent, effectively\naxiomatised and contains\nbasic arithmetic:", fontsize=9.3,
        color=INK, va="top")
ax.text(tx, 5.2, "First theorem (1931)", fontsize=10, color=PALETTE[4], weight="bold")
ax.text(tx, 4.95, "some sentence $G$ lies outside\nboth ovals: $T$ is incomplete.", fontsize=9, color=INK, va="top")
ax.text(tx, 3.75, "Second theorem (1931)", fontsize=10, color=PALETTE[3], weight="bold")
ax.text(tx, 3.5, "\"$T$ is consistent\" is one such\nsentence: $T$ cannot prove\nits own consistency.",
        fontsize=9, color=INK, va="top")
ax.text(tx, 1.95, "Adding $G$ as an axiom gives a\nnew theory $T'$ with its own $G'$.\nThe gap never closes.",
        fontsize=9, color=MUTED, va="top")
save(fig, "ch68-godel-map")

# ---- Goodstein sequences: digits of G_k(n) against step k


def bump(m, b):
    result, power = 0, 0
    while m:
        m, d = divmod(m, b)
        result += d * (b + 1) ** bump(power, b)
        power += 1
    return result


fig, ax = new_figure(figsize=(7.4, 4.2))
for n, col in zip((3, 4, 5, 6, 7), PALETTE):
    m, b, ys = n, 2, [math.log10(n)]
    for k in range(40):
        if m == 0:
            break
        m = bump(m, b) - 1
        b += 1
        ys.append(math.log10(m) if m > 0 else -0.3)
    ax.plot(range(len(ys)), ys, marker="o", ms=3, lw=1.6, color=col, label=f"start $n={n}$")
ax.set_xlabel("step $k$ (base is $k+2$)")
ax.set_ylabel(r"$\log_{10} G_k(n)$  (number of digits)")
ax.set_title("Goodstein sequences: every one reaches 0, eventually")
ax.annotate("$n=3$ dies at step 5", xy=(5, -0.3), xytext=(8, 1.2), fontsize=9, color=PALETTE[0],
            arrowprops=dict(arrowstyle="-|>", color=PALETTE[0], lw=1.2))
ax.annotate("$n=4$ climbs for about $3\\cdot 2^{402653211}$ steps\nbefore it finally falls to 0", xy=(32, 3.2),
            xytext=(39, 15), ha="right", fontsize=9, color=PALETTE[1],
            arrowprops=dict(arrowstyle="-|>", color=PALETTE[1], lw=1.2))
ax.legend(loc="upper left", ncol=1, fontsize=9)
ax.set_ylim(-0.6, 17)
save(fig, "ch68-goodstein")
