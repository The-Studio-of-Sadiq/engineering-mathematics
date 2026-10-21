"""Figures for Chapter 60 — Real Analysis: Making Calculus Rigorous."""
import sys, pathlib; sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
from style import plt, np, PALETTE, INK, MUTED, new_figure, diagram, node, box, arrow, save

BLUE, RED, GREEN, PURPLE, ORANGE = PALETTE[:5]

# ---------------------------------------------------------------- figure 1: epsilon-N picture
fig, ax = new_figure(figsize=(7.4, 4.0))
n = np.arange(1, 41)
a = 1 + (-1) ** n * 1.6 / n
L, eps = 1.0, 0.15
N = int(1.6 / eps)  # |a_n - 1| = 1.6/n < eps once n > 1.6/eps = 10.67
ax.axhspan(L - eps, L + eps, color=GREEN, alpha=0.15, lw=0)
ax.axhline(L, color=GREEN, lw=1.4, ls="--")
inside = n > N
ax.plot(n[~inside], a[~inside], "o", color=RED, ms=5, label=r"$n \leq N$: may stray outside")
ax.plot(n[inside], a[inside], "o", color=BLUE, ms=5, label=r"$n > N$: all inside the band")
ax.axvline(N + 0.5, color=INK, lw=1.2, ls=":")
ax.text(N + 0.9, 2.35, rf"$N = {N}$", color=INK, fontsize=11)
ax.annotate("", xy=(42, L + eps), xytext=(42, L - eps),
            arrowprops=dict(arrowstyle="<->", color=GREEN, lw=1.4))
ax.text(42.6, L, r"$2\varepsilon$", color=GREEN, va="center", fontsize=11)
ax.text(25, L + eps + 0.08, r"band $(L-\varepsilon,\ L+\varepsilon)$", color=GREEN, fontsize=10)
ax.set_xlabel("$n$"); ax.set_ylabel("$a_n$")
ax.set_xlim(0, 45); ax.set_ylim(-0.75, 2.6)
ax.set_title(r"$a_n = 1 + (-1)^n\,1.6/n \to 1$: with $\varepsilon = 0.15$, every $a_n$ with $n > 10$ is in the band", fontsize=11)
ax.legend(loc="lower right")
save(fig, "ch60-epsilon-n")

# ---------------------------------------------------------------- figure 2: Cauchy sequence in Q with a hole
fig, axs = new_figure(2, 1, figsize=(7.4, 3.9))
from fractions import Fraction as Fr
x = Fr(1); seq = []
for _ in range(5):
    seq.append(x); x = x / 2 + 1 / x
r2 = np.sqrt(2)

ax = axs[0]
ax.set_xlim(0.95, 1.55); ax.set_ylim(-0.6, 1.0); ax.axis("off")
ax.plot([0.95, 1.55], [0, 0], color=INK, lw=1.2)
for k, s in enumerate(seq[:2]):
    ax.plot(float(s), 0, "o", color=BLUE, ms=7, zorder=3)
    ax.text(float(s), 0.25, f"$x_{k}={s.numerator}/{s.denominator}$" if s.denominator > 1 else f"$x_{k}=1$",
            ha="center", color=INK, fontsize=10)
ax.plot(r2, 0, "o", mfc="white", mec=RED, ms=9, mew=2, zorder=4)
ax.text(r2, -0.45, r"$\sqrt{2}$: not in $\mathbb{Q}$", ha="center", color=RED, fontsize=10)
ax.text(r2 - 0.012, 0.25, r"$x_2, x_3,\ldots$ all in the box", ha="right", color=MUTED, fontsize=10)
ax.add_patch(plt.Rectangle((r2 - 0.006, -0.12), 0.012, 0.24, fill=False, ec=MUTED, lw=1, ls="--"))
ax.set_title(r"Newton's iterates $x_{k+1} = x_k/2 + 1/x_k$: every term rational", fontsize=11)

ax = axs[1]
lo, hi = r2 - 1.0e-3, r2 + 3.0e-3
ax.set_xlim(lo, hi); ax.set_ylim(-0.6, 1.0); ax.axis("off")
ax.plot([lo, hi], [0, 0], color=INK, lw=1.2)
ax.plot(float(seq[2]), 0, "o", color=BLUE, ms=7)
ax.text(float(seq[2]), 0.25, "$x_2=17/12$", ha="center", fontsize=10, color=INK)
ax.plot(float(seq[3]), 0, "o", color=BLUE, ms=7, zorder=3)
ax.text(float(seq[3]), 0.25, r"$x_3=577/408$, $x_4$, ... (within $2\times10^{-6}$)", ha="center", fontsize=10, color=INK)
ax.plot(r2, 0, "o", mfc="white", mec=RED, ms=9, mew=2, zorder=4)
ax.text(r2, -0.45, r"the gap they converge to", ha="center", color=RED, fontsize=10)
ax.set_title(r"zoom $\times 150$: the terms crowd together (Cauchy) but in $\mathbb{Q}$ there is nothing to land on",
             fontsize=11)
save(fig, "ch60-cauchy-hole")

# ---------------------------------------------------------------- figure 3: Weierstrass function with zooms
def W(x, a=0.5, b=5, terms=12):
    return sum(a ** k * np.cos(b ** k * np.pi * x) for k in range(terms))

fig, axs = new_figure(1, 3, figsize=(11.0, 3.5))
windows = [(-1.0, 1.0), (0.3, 0.3 + 0.4), (0.3, 0.3 + 0.08)]
cols = [BLUE, PURPLE, RED]
for i, ((x0, x1), ax) in enumerate(zip(windows, axs)):
    xs = np.linspace(x0, x1, 20001)
    ax.plot(xs, W(xs), color=cols[i], lw=0.7)
    ax.set_xlim(x0, x1)
    ax.set_title(["the whole graph", r"zoom $\times 5$", r"zoom $\times 25$"][i], fontsize=11)
    ax.set_xlabel("$x$")
    if i < 2:
        nx0, nx1 = windows[i + 1]
        ys = W(np.linspace(nx0, nx1, 4001))
        ax.add_patch(plt.Rectangle((nx0, ys.min()), nx1 - nx0, ys.max() - ys.min(),
                                   fill=False, ec=cols[i + 1], lw=1.4, zorder=4))
axs[0].set_ylabel(r"$W(x)=\sum_k 0.5^k\cos(5^k\pi x)$")
fig.suptitle("Continuous everywhere, differentiable nowhere: every zoom is as rough as the last",
             fontsize=12, fontweight="bold")
save(fig, "ch60-weierstrass")
