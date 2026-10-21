import sys, pathlib; sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
from style import plt, np, PALETTE, INK, MUTED, new_figure, diagram, node, box, arrow, save

BLUE, RED, GREEN, PURPLE, ORANGE, TEAL = PALETTE[:6]


def vec(ax, p, q, color, lw=2.2, ms=14):
    ax.annotate("", xy=q, xytext=p, zorder=5,
                arrowprops=dict(arrowstyle="-|>", color=color, lw=lw, mutation_scale=ms, shrinkA=0, shrinkB=0))


# ---- 1. Quotient space: cosets of a line U in the plane
fig, ax = new_figure(figsize=(7.6, 5.0))
xs = np.linspace(-4, 4, 2)
m = 0.5
for c in [-3, -2, -1, 1, 2, 3]:
    ax.plot(xs, m * xs + c, color=MUTED, lw=1.3, alpha=0.8)
ax.plot(xs, m * xs, color=BLUE, lw=3)
ax.text(-3.9, m * -3.9 - 0.45, "$U$ (the zero of $V/U$)", color=BLUE, fontsize=10, ha="left")
# two representatives of the same coset c = 2
v1 = np.array([-1.0, m * -1.0 + 2]); v2 = np.array([2.0, m * 2.0 + 2])
vec(ax, (0, 0), v1, RED); vec(ax, (0, 0), v2, RED)
ax.text(*(v1 + [-0.2, 0.25]), "$v$", color=RED, fontsize=12, ha="center")
ax.text(*(v2 + [0.05, 0.3]), "$v'$", color=RED, fontsize=12, ha="center")
ax.annotate("", xy=v2, xytext=v1, arrowprops=dict(arrowstyle="-|>", color=GREEN, lw=1.8, ls="--"))
ax.text(1.3, 2.15, "$v'-v\\in U$", color=GREEN, fontsize=10.5, ha="left")
ax.text(-3.85, 3.55, "grey lines: the cosets $v+U$", color=MUTED, fontsize=9.5, ha="left")
# a transversal W meeting every coset once
ax.plot([0, 0], [-4, 4], color=ORANGE, lw=1.6, ls=":")
ax.text(0.15, -3.75, "$W$ meets each coset once: $W\\cong V/U$", color=ORANGE, fontsize=9.5)
for c in [-3, -2, -1, 0, 1, 2, 3]:
    ax.plot(0, c, "o", color=ORANGE, ms=5, zorder=6)
ax.set_xlim(-4, 4); ax.set_ylim(-4, 4); ax.set_aspect("equal")
ax.set_xlabel("$x_1$"); ax.set_ylabel("$x_2$")
ax.set_title("$V=\\mathbb{R}^2$ sliced into cosets of $U$: each line is one point of $V/U$")
save(fig, "ch99-quotient")


# ---- 2. A basis and its dual basis, drawn as families of level lines
fig, ax = new_figure(figsize=(7.0, 5.4))
g = np.linspace(-3, 3, 2)
for k in range(-2, 4):          # phi1 = x - y = k
    ax.plot(g, g - k, color=PURPLE, lw=1.1, alpha=0.75)
for k in range(-2, 3):          # phi2 = y = k
    ax.plot(g, 0 * g + k, color=TEAL, lw=1.1, alpha=0.75)
for k in (0, 1):
    ax.plot(g, g - k, color=PURPLE, lw=2.4)
    ax.plot(g, 0 * g + k, color=TEAL, lw=2.4)
ax.text(-1.95, -1.9, r"$\varphi_1=0$", color=PURPLE, fontsize=10, ha="left")
ax.text(-0.75, -1.9, r"$\varphi_1=1$", color=PURPLE, fontsize=10, ha="left")
ax.text(-2.95, 0.07, r"$\varphi_2=0$", color=TEAL, fontsize=10, va="bottom")
ax.text(-2.95, 1.07, r"$\varphi_2=1$", color=TEAL, fontsize=10, va="bottom")
vec(ax, (0, 0), (1, 0), RED, lw=2.6); vec(ax, (0, 0), (1, 1), BLUE, lw=2.6)
ax.text(1.12, -0.32, r"$e_1=(1,0)$", color=RED, fontsize=11)
ax.text(0.55, 1.25, r"$e_2=(1,1)$", color=BLUE, fontsize=11, ha="right")
ax.set_xlim(-3, 3); ax.set_ylim(-2.4, 2.6); ax.set_aspect("equal")
ax.set_title(r"dual basis: $\varphi_1=x-y$, $\varphi_2=y$, so $\varphi_i(e_j)=\delta_{ij}$")
ax.set_xlabel(r"a covector is a stack of level lines; $\varphi(v)$ counts the lines $v$ crosses",
              fontsize=9.5, color=MUTED)
save(fig, "ch99-dual-basis")


# ---- 3. Jordan chains of a nilpotent map with blocks of sizes 3, 2, 1
fig, ax = diagram(figsize=(9.0, 4.8))
ax.set_xlim(0, 12); ax.set_ylim(0, 6.4)
cols = [(2.0, 3), (4.4, 2), (6.8, 1)]
ylev = {3: 5.0, 2: 3.5, 1: 2.0}
for x, size in cols:
    for lvl in range(size, 0, -1):
        lab = "$v$" if lvl == size else ("$N^{%d}v$" % (size - lvl) if size - lvl > 1 else "$Nv$")
        node(ax, x, ylev[lvl], lab, color=[BLUE, GREEN, ORANGE][[3, 2, 1].index(size)], r=0.42, fontsize=9.5)
    for lvl in range(size, 1, -1):
        arrow(ax, x, ylev[lvl], x, ylev[lvl - 1], color=INK, text=None, shrink=25)
    arrow(ax, x, ylev[1], x, 0.55, color=MUTED, shrink=25)
    ax.text(x, 0.75, "0", ha="center", va="center", fontsize=11, color=MUTED)
    ax.text(x, 0.15, f"block of size {size}", ha="center", fontsize=9, color=MUTED)
for lvl, txt in [(1, r"$\ker N$: dim 3"), (2, r"$\ker N^2$: dim 5"), (3, r"$\ker N^3$: dim 6")]:
    ax.plot([1.2, 7.6], [ylev[lvl] - 0.62] * 2, color="#d1d5db", lw=0.8, ls="--")
    ax.text(8.0, ylev[lvl], txt, fontsize=10, va="center", color=INK)
ax.text(8.0, 5.9, r"every arrow is $N$; each column is one chain", fontsize=9.5, color=MUTED)
ax.text(8.0, 0.75, r"blocks of size $\geq k$ $=\ \dim\ker N^k-\dim\ker N^{k-1}$", fontsize=9.5, color=PURPLE)
save(fig, "ch99-jordan-chain")
