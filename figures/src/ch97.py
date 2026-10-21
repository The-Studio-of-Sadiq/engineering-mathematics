"""Figures for Chapter 97 — Measure Theory and Modern Probability."""
import sys, pathlib; sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
from style import plt, np, PALETTE, INK, MUTED, new_figure, diagram, node, box, arrow, save

BLUE, RED, GREEN, PURPLE, ORANGE, TEAL, GREY = PALETTE

f = lambda x: 0.55 + 0.28 * np.sin(2 * np.pi * x) + 0.12 * np.sin(6 * np.pi * x + 0.6)
xx = np.linspace(0, 1, 1200)

# ---------------------------------------------------------------- figure 1: Riemann vs Lebesgue
fig, axs = new_figure(1, 2, figsize=(9.4, 3.9), sharey=True)
ax = axs[0]
n = 12
edges = np.linspace(0, 1, n + 1)
for k in range(n):
    xm = 0.5 * (edges[k] + edges[k + 1])
    ax.add_patch(plt.Rectangle((edges[k], 0), 1 / n, f(xm), fc=BLUE, alpha=0.22, ec=BLUE, lw=0.8))
ax.plot(xx, f(xx), color=INK, lw=2)
ax.set_title("Riemann: slice the domain\n(sum of height × width)", fontsize=10.5)
ax.set_xlabel("$x$"); ax.set_ylabel("$f(x)$")

ax = axs[1]
levels = np.linspace(0, 0.95, 20)
dy = levels[1] - levels[0]
fx = f(xx)
for k, y in enumerate(levels):
    above = fx > y + dy / 2
    meas = above.mean()
    ax.add_patch(plt.Rectangle((0, y), meas, dy, fc=GREEN, alpha=0.20, ec=GREEN, lw=0.6))
ax.plot(xx, fx, color=INK, lw=2)
y0 = 0.62
ax.axhline(y0, color=RED, lw=1.4, ls="--")
mask = fx > y0
# mark the set {f > y0} on the x-axis
starts = np.where(np.diff(np.r_[0, mask.astype(int), 0]) == 1)[0]
ends = np.where(np.diff(np.r_[0, mask.astype(int), 0]) == -1)[0] - 1
for s_, e_ in zip(starts, ends):
    ax.plot([xx[s_], xx[e_]], [0.012, 0.012], color=RED, lw=6, solid_capstyle="butt")
ax.text(0.98, y0 + 0.03, "level $y$", color=RED, ha="right", fontsize=9.5)
ax.text(0.5, 0.07, r"the set $\{f>y\}$: its measure $m\{f>y\}$ is the slice's length", color=RED,
        ha="center", fontsize=9)
ax.set_title("Lebesgue: slice the range\n" r"$\int f\,dm=\int_0^\infty m\{f>y\}\,dy$", fontsize=10.5)
ax.set_xlabel("$x$  (green bars: length of $\\{f>y\\}$ at each level, stacked)")
for a in axs:
    a.set_xlim(0, 1); a.set_ylim(0, 1.0)
save(fig, "ch97-riemann-vs-lebesgue")


# ---------------------------------------------------------------- figure 2: Cantor set and Cantor function
def cantor_intervals(level):
    iv = [(0.0, 1.0)]
    for _ in range(level):
        iv = [p for a, b in iv for p in ((a, a + (b - a) / 3), (b - (b - a) / 3, b))]
    return iv


def cantor_function(x, depth=30):
    x = np.asarray(x, float).copy(); out = np.zeros_like(x); scale = 0.5
    done = np.zeros_like(x, bool)
    for _ in range(depth):
        d = np.floor(x * 3).clip(0, 2); x = x * 3 - d
        out += np.where(~done & (d == 2), scale, 0.0)
        hit = ~done & (d == 1)
        out += np.where(hit, scale, 0.0); done |= hit
        scale /= 2
    return out


fig, axs = new_figure(1, 2, figsize=(9.6, 3.9), gridspec_kw={"width_ratios": [1.05, 1]})
ax = axs[0]
for lev in range(6):
    yv = 5 - lev
    for a_, b_ in cantor_intervals(lev):
        ax.plot([a_, b_], [yv, yv], color=PURPLE, lw=7, solid_capstyle="butt")
    ax.text(1.04, yv, f"$C_{lev}$: total length $(2/3)^{lev}\\approx{(2/3)**lev:.3f}$", va="center", fontsize=9)
ax.set_xlim(0, 1.75); ax.set_ylim(-0.6, 5.6); ax.axis("off")
ax.set_title("Remove middle thirds forever:\nlength $\\to 0$, yet uncountably many points remain", fontsize=10.5)

ax = axs[1]
x = np.linspace(0, 1, 6000)
ax.plot(x, cantor_function(x), color=PURPLE, lw=1.8)
ax.axvspan(1 / 3, 2 / 3, color=GREEN, alpha=0.10, lw=0)
ax.text(0.5, 0.33, "flat on every\nremoved interval", ha="center", fontsize=9, color=GREEN)
ax.text(0.03, 0.88, "climbs from 0 to 1 using\nonly the Cantor set,\nwhich has length 0", fontsize=9, color=INK)
ax.set_xlabel("$x$"); ax.set_ylabel("$F(x)$")
ax.set_xlim(0, 1); ax.set_ylim(0, 1.02)
ax.set_title("The Cantor function: a continuous CDF\nwith no density and no jumps", fontsize=10.5)
save(fig, "ch97-cantor")
