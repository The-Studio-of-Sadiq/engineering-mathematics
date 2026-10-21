import sys, pathlib; sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
from style import plt, np, PALETTE, INK, MUTED, new_figure, diagram, node, box, arrow, save
from matplotlib.patches import FancyBboxPatch

BLUE, RED, GREEN, PURPLE, ORANGE = PALETTE[0], PALETTE[1], PALETTE[2], PALETTE[3], PALETTE[4]


def harrow(ax, p, q, text=None, color=INK, curve=0.0, off=(0, 0.2), shrink=20, ls="-", lw=1.6, fs=10,
           tcolor=None, shrinkB=None):
    """Arrow from p to q with a label placed at the midpoint plus an offset."""
    ax.annotate("", xy=q, xytext=p, zorder=2,
                arrowprops=dict(arrowstyle="-|>", color=color, lw=lw, shrinkA=shrink, shrinkB=shrink if shrinkB is None else shrinkB,
                                linestyle=ls, connectionstyle=f"arc3,rad={curve}"))
    if text:
        mx, my = (p[0] + q[0]) / 2, (p[1] + q[1]) / 2
        ax.text(mx + off[0], my + off[1], text, ha="center", va="center", fontsize=fs,
                color=tcolor or color, zorder=4)


def region(ax, x, y, w, h, title, color):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.02,rounding_size=0.35",
                                fc=color + "10", ec=color, lw=1.6, zorder=0))
    ax.text(x + w / 2, y + h - 0.3, title, ha="center", va="center", fontsize=11, color=color,
            weight="bold")


# ---- 1. A functor F : C -> D
fig, ax = diagram(figsize=(9.0, 4.6))
ax.set_xlim(0, 12); ax.set_ylim(0, 6.2)
region(ax, 0.2, 0.5, 4.6, 5.0, r"category $\mathcal{C}$", BLUE)
region(ax, 7.2, 0.5, 4.6, 5.0, r"category $\mathcal{D}$", GREEN)
C = {"A": (1.2, 1.6), "B": (2.6, 4.4), "C": (4.0, 1.6)}
D = {"A": (8.0, 1.6), "B": (9.4, 4.4), "C": (10.8, 1.6)}
for k in C:
    node(ax, *C[k], f"${k}$", color=BLUE, r=0.38, fontsize=11)
    node(ax, *D[k], f"$F{k}$", color=GREEN, r=0.45, fontsize=10)
harrow(ax, C["A"], C["B"], "$f$", color=BLUE, off=(-0.35, 0.05), shrink=24)
harrow(ax, C["B"], C["C"], "$g$", color=BLUE, off=(0.35, 0.05), shrink=24)
harrow(ax, C["A"], C["C"], r"$g\circ f$", color=BLUE, off=(0, -0.3), shrink=24)
harrow(ax, D["A"], D["B"], "$Ff$", color=GREEN, off=(-0.45, 0.05), shrink=28)
harrow(ax, D["B"], D["C"], "$Fg$", color=GREEN, off=(0.45, 0.05), shrink=28)
harrow(ax, D["A"], D["C"], r"$F(g\circ f)=Fg\circ Ff$", color=GREEN, off=(0, -0.32), shrink=28)
for k, curve in [("A", 0.3), ("B", -0.12), ("C", -0.2)]:
    harrow(ax, C[k], D[k], None, color=MUTED, ls="--", curve=curve, lw=1.2, shrink=24, shrinkB=28)
ax.text(6.0, 5.95, r"$F$ sends objects to objects", ha="center", fontsize=9.5, color=MUTED)
ax.text(6.0, 5.6, "and arrows to arrows", ha="center", fontsize=9.5, color=MUTED)
ax.text(6.0, 0.15, r"laws: $F(\mathrm{id}_A)=\mathrm{id}_{FA}$ and $F(g\circ f)=Fg\circ Ff$;"
        r"   example: $F=$ List, $Ff=$ map $f$", ha="center", fontsize=9.5, color=INK)
save(fig, "ch62-functor")


# ---- 2. A naturality square, abstract and concrete
fig, axes = plt.subplots(1, 2, figsize=(10.0, 4.2), constrained_layout=True)
for ax in axes:
    ax.set_xlim(0, 6); ax.set_ylim(0, 5); ax.set_aspect("equal"); ax.axis("off")
ax = axes[0]
P = {"FA": (1.0, 3.8), "FB": (5.0, 3.8), "GA": (1.0, 1.0), "GB": (5.0, 1.0)}
lab = {"FA": "$FA$", "FB": "$FB$", "GA": "$GA$", "GB": "$GB$"}
for k, (x, y) in P.items():
    ax.text(x, y, lab[k], ha="center", va="center", fontsize=13, color=INK)
harrow(ax, P["FA"], P["FB"], "$Ff$", color=GREEN, shrink=18)
harrow(ax, P["GA"], P["GB"], "$Gf$", color=PURPLE, off=(0, -0.28), shrink=18)
harrow(ax, P["FA"], P["GA"], r"$\eta_A$", color=RED, off=(-0.4, 0), shrink=16)
harrow(ax, P["FB"], P["GB"], r"$\eta_B$", color=RED, off=(0.4, 0), shrink=16)
ax.text(3.0, 2.4, "commutes:", ha="center", fontsize=9.5, color=MUTED)
ax.text(3.0, 2.0, r"$\eta_B\circ Ff = Gf\circ\eta_A$", ha="center", fontsize=10.5, color=INK)
ax.set_title(r"natural transformation $\eta: F\Rightarrow G$", fontsize=11)

ax = axes[1]
Q = {"FA": (1.0, 3.8), "FB": (5.0, 3.8), "GA": (1.0, 1.0), "GB": (5.0, 1.0)}
txt = {"FA": "[1, 2, 4]", "FB": "[F, T, T]", "GA": "[4, 2, 1]", "GB": "[T, T, F]"}
for k, (x, y) in Q.items():
    ax.add_patch(FancyBboxPatch((x - 0.8, y - 0.3), 1.6, 0.6, boxstyle="round,pad=0.02,rounding_size=0.1",
                                fc="#f9fafb", ec=MUTED, lw=1.0))
    ax.text(x, y, txt[k], ha="center", va="center", fontsize=10.5, family="monospace", color=INK)
harrow(ax, Q["FA"], Q["FB"], "map even", color=GREEN, shrink=30, fs=9.5)
harrow(ax, Q["GA"], Q["GB"], "map even", color=PURPLE, off=(0, -0.3), shrink=30, fs=9.5)
harrow(ax, Q["FA"], Q["GA"], "reverse", color=RED, off=(-0.55, 0), shrink=14, fs=9.5)
harrow(ax, Q["FB"], Q["GB"], "reverse", color=RED, off=(0.55, 0), shrink=14, fs=9.5)
ax.text(3.0, 2.4, "both paths give", ha="center", fontsize=9.5, color=MUTED)
ax.text(3.0, 2.0, "[T, T, F]", ha="center", fontsize=10.5, family="monospace", color=INK)
ax.set_title(r"instance: $F=G=$ List, $\eta=$ reverse, $f=$ even", fontsize=11)
save(fig, "ch62-naturality")


# ---- 3. The universal property of a product, with four instances
fig, ax = diagram(figsize=(9.6, 5.4))
ax.set_xlim(0, 12); ax.set_ylim(0, 6.8)
X, PR, A, B = (3.2, 6.0), (3.2, 3.4), (0.9, 1.4), (5.5, 1.4)
node(ax, *X, "$X$", color=ORANGE, r=0.42, fontsize=12)
box(ax, *PR, r"$A\times B$", w=1.5, h=0.7, color=BLUE, fontsize=12)
node(ax, *A, "$A$", color=BLUE, r=0.42, fontsize=12)
node(ax, *B, "$B$", color=BLUE, r=0.42, fontsize=12)
harrow(ax, X, A, "$f$", color=ORANGE, curve=0.25, off=(-0.85, 0.2), shrink=27)
harrow(ax, X, B, "$g$", color=ORANGE, curve=-0.25, off=(0.85, 0.2), shrink=27)
harrow(ax, X, PR, r"$\exists!\ \langle f,g\rangle$", color=RED, ls="--", off=(0.95, 0.0), shrink=27, shrinkB=22)
harrow(ax, PR, A, r"$\pi_1$", color=BLUE, off=(-0.42, 0.35), shrink=27)
harrow(ax, PR, B, r"$\pi_2$", color=BLUE, off=(0.42, 0.35), shrink=27)
ax.text(3.2, 0.35, r"$\pi_1\circ\langle f,g\rangle=f,\quad \pi_2\circ\langle f,g\rangle=g$",
        ha="center", fontsize=10.5, color=INK)
rows = [("category", r"product $A\times B$", BLUE),
        ("Sets", "pairs (a, b)", INK),
        ("Vector spaces", r"direct sum $A\oplus B$", INK),
        ("Divisibility of integers", r"$\gcd(a,b)$", INK),
        ("Real numbers with $\\leq$", r"$\min(a,b)$", INK),
        ("Logic (implication)", r"$a$ AND $b$", INK),
        ("Program types", "tuple (A, B)", INK)]
for i, (c, p, col) in enumerate(rows):
    y = 6.0 - 0.72 * i
    ax.text(7.2, y, c, ha="left", va="center", fontsize=10, color=col, weight="bold" if i == 0 else None)
    ax.text(10.4, y, p, ha="left", va="center", fontsize=10, color=col, weight="bold" if i == 0 else None)
ax.plot([7.1, 12.3], [5.62, 5.62], color=MUTED, lw=0.8)
ax.text(9.7, 0.65, "one definition, six familiar constructions", ha="center", fontsize=9.5, color=MUTED)
save(fig, "ch62-product")
