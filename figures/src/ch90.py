"""Figures for Chapter 90 — Spectral Theory: Understanding Operators Through Their Modes."""
import sys, pathlib; sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
from style import plt, np, PALETTE, INK, MUTED, new_figure, diagram, node, box, arrow, save

BLUE, RED, GREEN, PURPLE, ORANGE, TEAL, GREY = PALETTE
pi = np.pi

# ---------------------------------------------------------------- figure 1: modes of a square drum
g = np.linspace(0, 1, 201)
X, Y = np.meshgrid(g, g)


def mode(m, n):
    return np.sin(m * pi * X) * np.sin(n * pi * Y)


panels = [
    ("(1,1)", mode(1, 1), 2), ("(2,1)", mode(2, 1), 5), ("(1,2)", mode(1, 2), 5), ("(2,2)", mode(2, 2), 8),
    ("(2,1)+(1,2)", mode(2, 1) + mode(1, 2), 5), ("(2,1)$-$(1,2)", mode(2, 1) - mode(1, 2), 5),
    ("(3,1)+(1,3)", mode(3, 1) + mode(1, 3), 10), ("(3,2)", mode(3, 2), 13),
]
fig, axs = plt.subplots(2, 4, figsize=(10.4, 5.6), constrained_layout=True)
for ax, (name, Z, lam) in zip(axs.ravel(), panels):
    Z = Z / np.abs(Z).max()
    ax.pcolormesh(X, Y, Z, cmap="RdBu_r", vmin=-1, vmax=1, shading="auto")
    ax.contour(X, Y, Z, levels=[0], colors=INK, linewidths=1.2)
    ax.set_aspect("equal"); ax.set_xticks([]); ax.set_yticks([]); ax.grid(False)
    for s in ax.spines.values():
        s.set_visible(True); s.set_color(INK)
    ax.set_title(f"{name}\n" + rf"$\lambda={lam}\pi^2$", fontsize=10)
fig.suptitle("Modes of a square drum: red and blue move in opposite directions, black lines stay still",
             fontsize=11, fontweight="bold")
save(fig, "ch90-drum-modes")

# ---------------------------------------------------------------- figure 2: counting function vs Weyl's law
fig, axs = new_figure(1, 2, figsize=(10.6, 4.2))
ax = axs[0]
ev = np.sort(np.array([pi**2 * (m * m + n * n) for m in range(1, 60) for n in range(1, 60)]))
lam = np.linspace(0, 1500, 3001)
Ncont = np.searchsorted(ev, lam, side="right")
ax.step(lam, Ncont, where="post", color=BLUE, lw=1.8, label=r"exact count $N(\lambda)$, unit square")
n = 10
h = 1 / (n + 1)
s = 4 / h**2 * np.sin(np.arange(1, n + 1) * pi * h / 2) ** 2
evh = np.sort((s[:, None] + s[None, :]).ravel())
ax.step(lam, np.searchsorted(evh, lam, side="right"), where="post", color=ORANGE, lw=1.5,
        label=r"discrete Laplacian, $10\times10$ grid")
ax.plot(lam, lam / (4 * pi), color=GREY, ls="--", lw=1.6, label=r"Weyl: $\frac{A}{4\pi}\lambda$")
ax.plot(lam, lam / (4 * pi) - np.sqrt(lam) / pi, color=RED, lw=1.6,
        label=r"two terms: $\frac{A}{4\pi}\lambda-\frac{P}{4\pi}\sqrt{\lambda}$")
ax.set_xlim(0, 1500); ax.set_ylim(0, 130)
ax.set_xlabel(r"$\lambda$"); ax.set_ylabel(r"number of eigenvalues $\leq\lambda$")
ax.legend(fontsize=8.6, loc="upper left")
ax.set_title("Low frequencies: the staircase hugs Weyl's law", fontsize=11)

ax = axs[1]
for n, col in [(10, ORANGE), (20, GREEN), (40, PURPLE)]:
    h = 1 / (n + 1)
    s = 4 / h**2 * np.sin(np.arange(1, n + 1) * pi * h / 2) ** 2
    evh = np.sort((s[:, None] + s[None, :]).ravel())
    ax.plot(evh * h**2, np.arange(1, n * n + 1) / n**2, color=col, lw=1.8, label=f"grid ${n}\\times{n}$")
xx = np.linspace(0, 8, 400)
ax.plot(xx, xx / (4 * pi), color=GREY, ls="--", lw=1.6, label=r"Weyl's law, $\lambda h^2/(4\pi)$")
ax.axhline(1, color=MUTED, lw=0.8, ls=":")
ax.set_xlim(0, 8.2); ax.set_ylim(0, 1.08)
ax.set_xlabel(r"scaled eigenvalue $\lambda h^2$"); ax.set_ylabel(r"fraction of all $n^2$ eigenvalues")
ax.text(5.3, 0.12, "the grid runs out of\nmodes at $\\lambda=8/h^2$", fontsize=9.5, color=MUTED)
ax.legend(fontsize=8.6, loc="upper left")
ax.set_title("High frequencies: every grid departs from Weyl", fontsize=11)
save(fig, "ch90-weyl-law")

# ---------------------------------------------------------------- figure 3: three kinds of spectrum
fig, ax = plt.subplots(figsize=(9.0, 3.6), constrained_layout=True)
ax.axis("off")
ax.set_xlim(-1.35, 1.25); ax.set_ylim(-0.3, 3.2)
rows = [
    (2.6, "string on $[0,1]$", r"$-\frac{d^2}{dx^2}$, $u(0)=u(1)=0$", "pure point: $k^2\\pi^2$, running off to $\\infty$"),
    (1.5, "free particle on $\\mathbb{R}$", r"$-\frac{d^2}{dx^2}$ on $L^2(\mathbb{R})$", "purely continuous: all of $[0,\\infty)$, no eigenvectors"),
    (0.4, "hydrogen atom", r"$-\frac{1}{2}\nabla^2-\frac{1}{r}$", "both: $-1/(2n^2)$ piling up at $0$, then a continuum"),
]
for y, name, op, desc in rows:
    ax.plot([-0.55, 1.15], [y, y], color=INK, lw=1)
    ax.text(-1.33, y + 0.12, name, fontsize=10.5, fontweight="bold", color=INK, va="bottom")
    ax.text(-1.33, y - 0.08, op, fontsize=10, color=INK, va="top")
    ax.text(0.3, y - 0.22, desc, fontsize=9.5, color=MUTED, ha="center", va="top")
# string: points at k^2 pi^2 scaled
y = 2.6
for k in range(1, 4):
    xv = -0.5 + 0.18 * k**2
    ax.plot(xv, y, "o", color=BLUE, ms=8, zorder=3)
ax.text(1.17, y, r"$\cdots$", fontsize=12, va="center", color=BLUE)
ax.plot(-0.5, y, "|", color=INK, ms=14)
ax.text(-0.5, y + 0.12, "0", ha="center", fontsize=9, color=INK)
# free particle
y = 1.5
ax.plot([-0.5, 1.15], [y, y], color=RED, lw=6, solid_capstyle="butt", alpha=0.8)
ax.plot(-0.5, y, "|", color=INK, ms=14)
ax.text(-0.5, y + 0.12, "0", ha="center", fontsize=9, color=INK)
# hydrogen
y = 0.4
x0 = 0.25
for nq in range(1, 30):
    xv = x0 - 0.7 / nq**2
    ax.plot(xv, y, "o", color=BLUE, ms=max(2.5, 8 - nq), zorder=3)
ax.plot([x0, 1.15], [y, y], color=RED, lw=6, solid_capstyle="butt", alpha=0.8)
ax.plot(x0, y, "|", color=INK, ms=14)
ax.text(x0, y + 0.12, "0", ha="center", fontsize=9, color=INK)
ax.plot([], [], "o", color=BLUE, label="eigenvalues (point spectrum)")
ax.plot([], [], color=RED, lw=6, alpha=0.8, label="continuous spectrum")
ax.legend(loc="upper right", fontsize=9, ncol=2, bbox_to_anchor=(1.0, 1.08))
save(fig, "ch90-spectrum-types")
