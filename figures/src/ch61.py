"""Figures for Chapter 61 — Functional Analysis: Infinite-Dimensional Thinking."""
import sys, pathlib; sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
from style import plt, np, PALETTE, INK, MUTED, new_figure, diagram, node, box, arrow, save

BLUE, RED, GREEN, PURPLE, ORANGE = PALETTE[:5]

# ---------------------------------------------------------------- figure 1: unit balls of p-norms
fig, ax = new_figure(figsize=(5.6, 5.4))
t = np.linspace(0, 2 * np.pi, 2001)
c, s = np.cos(t), np.sin(t)
for p, col, lab in [(1, RED, r"$p=1$ (diamond)"), (2, BLUE, r"$p=2$ (circle)"), (4, GREEN, r"$p=4$")]:
    r = (np.abs(c) ** p + np.abs(s) ** p) ** (-1 / p)
    ax.plot(r * c, r * s, color=col, label=lab)
ax.plot([1, 1, -1, -1, 1], [1, -1, -1, 1, 1], color=PURPLE, label=r"$p=\infty$ (square)")
p = 0.5
r = (np.abs(c) ** p + np.abs(s) ** p) ** (-1 / p)
ax.plot(r * c, r * s, color=MUTED, ls="--", lw=1.6, label=r"$p=\frac{1}{2}$: not convex, not a norm")
ax.axhline(0, color=INK, lw=0.8); ax.axvline(0, color=INK, lw=0.8)
ax.set_aspect("equal"); ax.set_xlim(-1.35, 1.35); ax.set_ylim(-1.35, 1.65)
ax.set_xlabel("$x_1$"); ax.set_ylabel("$x_2$")
ax.legend(loc="upper center", ncol=2, fontsize=9, bbox_to_anchor=(0.5, 1.0))
ax.set_title(r"Unit balls $\Vert\mathbf{x}\Vert_p \leq 1$ in $\mathbb{R}^2$")
save(fig, "ch61-unit-balls")

# ---------------------------------------------------------------- figure 2: L2 convergence without uniform convergence
fig, axs = new_figure(1, 2, figsize=(10.0, 3.8))
ax = axs[0]
x = np.linspace(0, 1, 1001)
ns = [1, 2, 5, 10, 30, 100]
cmap = plt.get_cmap("viridis")
for i, n in enumerate(ns):
    ax.plot(x, x ** n, color=cmap(i / (len(ns) - 1) * 0.9), lw=1.8, label=f"$n={n}$")
ax.plot([0, 1], [0, 0], color=RED, lw=2.4, ls="--", zorder=1)
ax.set_ylim(-0.05, 1.05)
ax.plot(1, 1, "o", color=RED, ms=6); ax.plot(1, 0, "o", mfc="white", mec=RED, ms=6)
ax.text(0.03, 0.09, "pointwise limit: 0 on [0,1), 1 at x=1", color=RED, fontsize=9.5, va="bottom")
ax.set_xlabel("$x$"); ax.set_ylabel("$f_n(x) = x^n$")
ax.set_title(r"$f_n(x)=x^n$ on $[0,1]$", fontsize=11)
ax.legend(fontsize=8.5, loc="upper left", bbox_to_anchor=(0.0, 0.92))

ax = axs[1]
n = np.arange(1, 201)
ax.plot(n, np.ones_like(n, dtype=float), color=PURPLE, label=r"uniform distance $\sup_x | f_n - 0| = 1$")
ax.plot(n, 1 / np.sqrt(2 * n + 1), color=BLUE, label=r"$L^2$ distance $\Vert f_n\Vert_2 = 1/\sqrt{2n+1}$")
ax.set_xscale("log"); ax.set_yscale("log")
ax.set_ylim(0.03, 8)
ax.set_xlabel("$n$"); ax.set_ylabel("distance from the zero function")
ax.set_title(r"converges to $0$ in $L^2$, never uniformly", fontsize=11)
ax.legend(fontsize=9, loc="upper right")
save(fig, "ch61-l2-vs-uniform")

# ---------------------------------------------------------------- figure 3: hierarchy of spaces
fig, ax = diagram(figsize=(8.4, 4.4))
from matplotlib.patches import FancyBboxPatch
levels = [
    (0.2, 0.15, 9.6, 4.9, GREEN, "vector spaces", "polynomials of any degree, all sequences"),
    (0.8, 0.4, 8.4, 3.75, ORANGE, "normed spaces", r"$C[0,1]$ with $\Vert f\Vert_1=\int| f|$ (has holes)"),
    (1.4, 0.65, 7.2, 2.6, BLUE, "Banach spaces: normed + complete", r"$C[0,1]$ with sup norm, $L^1$, $\ell^\infty$"),
    (2.0, 0.9, 6.0, 1.45, PURPLE, "Hilbert spaces: inner product + complete", r"$\mathbb{R}^n$, $\ell^2$, $L^2$: the home of Fourier and quantum states"),
]
for (x0, y0, w, h, col, name, ex) in levels:
    ax.add_patch(FancyBboxPatch((x0, y0), w, h, boxstyle="round,pad=0.02,rounding_size=0.25",
                                fc=col, alpha=0.08, ec=col, lw=1.8))
    ax.add_patch(FancyBboxPatch((x0, y0), w, h, boxstyle="round,pad=0.02,rounding_size=0.25",
                                fc="none", ec=col, lw=1.8))
    ax.text(x0 + w / 2, y0 + h - 0.32, name, ha="center", va="center", fontsize=10.5, fontweight="bold", color=col)
    ax.text(x0 + w / 2, y0 + h - 0.68, ex, ha="center", va="center", fontsize=9, color=INK)
ax.set_ylim(0, 5.2)
save(fig, "ch61-space-hierarchy")
