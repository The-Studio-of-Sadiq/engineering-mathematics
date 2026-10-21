"""Figures for Chapter 55 — Optimal Transport and Wasserstein Geometry."""
import sys, pathlib; sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
from style import plt, np, PALETTE, INK, MUTED, new_figure, diagram, node, box, arrow, save
from scipy.stats import norm

BLUE, RED, GREEN, PURPLE, ORANGE, TEAL, GREY = PALETTE

# ------------------------------------------------- figure 1: the monotone transport map in 1-D
x = np.linspace(-5, 6, 2000)
src = 0.6 * norm.pdf(x, -2.0, 0.6) + 0.4 * norm.pdf(x, 1.0, 0.8)     # source: bimodal
tgt = norm.pdf(x, 2.0, 1.1)                                           # target: one bump
F = np.cumsum(src); F /= F[-1]
G = np.cumsum(tgt); G /= G[-1]
us = np.linspace(0.04, 0.96, 15)
xq = np.interp(us, F, x)
yq = np.interp(us, G, x)

fig, ax = new_figure(figsize=(8.6, 4.6))
ax.grid(False)
top, bot, s = 2.2, 0.0, 1.6
ax.fill_between(x, top, top + s * src, color=BLUE, alpha=0.25)
ax.plot(x, top + s * src, color=BLUE, lw=2)
ax.fill_between(x, bot, bot - s * tgt, color=ORANGE, alpha=0.25)
ax.plot(x, bot - s * tgt, color=ORANGE, lw=2)
ax.axhline(top, color=INK, lw=0.8); ax.axhline(bot, color=INK, lw=0.8)
for a, b_ in zip(xq, yq):
    ax.annotate("", xy=(b_, bot + 0.04), xytext=(a, top - 0.04),
                arrowprops=dict(arrowstyle="-|>", color=MUTED, lw=1.1, shrinkA=0, shrinkB=0))
    ax.plot(a, top, "o", color=BLUE, ms=3.5, zorder=5)
    ax.plot(b_, bot, "o", color=ORANGE, ms=3.5, zorder=5)
ax.text(-4.9, top + 0.55, r"source $\mu$" + "\n(two bumps)", color=BLUE, fontsize=10, va="center")
ax.text(-4.9, bot - 0.35, r"target $\nu$", color=ORANGE, fontsize=10, va="center")
ax.text(5.9, 1.1, r"$T=G^{-1}\circ F$:" + "\nthe $k$-th quantile\ngoes to the $k$-th\nquantile; no\narrows cross",
        ha="right", va="center", fontsize=9.5, color=INK)
ax.set_xlim(-5, 6); ax.set_ylim(-0.75, 3.0)
ax.set_yticks([]); ax.spines["left"].set_visible(False)
ax.set_xlabel("$x$")
save(fig, "ch55-transport-map")

# ------------------------------------------------- figure 2: Wasserstein vs KL, TV, JS
d = np.linspace(0, 4, 400)
fig, axs = new_figure(1, 2, figsize=(10.5, 4.2))
ax = axs[0]
ax.plot(d, d, color=BLUE, label=r"$W_1=W_2=d$")
ax.plot(d, d ** 2 / 2, color=RED, label=r"KL $=d^2/2$")
ax.plot(d, 2 * norm.cdf(d / 2) - 1, color=GREEN, label=r"TV $=2\Phi(d/2)-1$")
# JS divergence between N(0,1) and N(d,1), numerically (nats)
z = np.linspace(-12, 16, 6000); dz = z[1] - z[0]
js = []
for di in d:
    p, q = norm.pdf(z), norm.pdf(z, di)
    m = (p + q) / 2
    with np.errstate(divide="ignore", invalid="ignore"):
        t = 0.5 * np.where(p > 0, p * np.log(p / m), 0) + 0.5 * np.where(q > 0, q * np.log(q / m), 0)
    js.append(np.nansum(t) * dz)
ax.plot(d, js, color=PURPLE, ls="--", label=r"JS (saturates at $\ln 2$)")
ax.set_xlim(0, 4); ax.set_ylim(0, 4.2)
ax.set_xlabel("shift $d$"); ax.set_ylabel("discrepancy")
ax.set_title(r"Gaussians $\mathcal{N}(0,1)$ vs $\mathcal{N}(d,1)$", fontsize=11)
ax.legend(fontsize=9, loc="upper left")

ax = axs[1]
ax.plot(d, d, color=BLUE, label=r"$W_1=d$")
ax.plot(d, np.minimum(d, 1), color=GREEN, label=r"TV $=\min(d,1)$")
ax.plot(d, np.log(2) * np.minimum(d, 1), color=PURPLE, ls="--", label=r"JS $=\ln 2\cdot\min(d,1)$")
ax.annotate(r"KL $=\infty$ for every $d>0$", xy=(0.05, 4.05), xytext=(0.6, 3.3), color=RED, fontsize=10,
            arrowprops=dict(arrowstyle="-|>", color=RED, lw=1))
ax.axvspan(1, 4, color=GREY, alpha=0.06)
ax.text(2.6, 1.35, "disjoint supports:\nTV and JS are flat,\nso their gradient is zero", ha="center",
        fontsize=9, color=MUTED)
ax.set_xlim(0, 4); ax.set_ylim(0, 4.2)
ax.set_xlabel("shift $d$")
ax.set_title(r"uniform boxes $[0,1]$ vs $[d,d+1]$", fontsize=11)
ax.legend(fontsize=9, loc="upper right")
save(fig, "ch55-divergences")

# ------------------------------------------------- figure 3: Sinkhorn couplings for several epsilon
n = 60
grid = np.linspace(0, 1, n)
a = 0.5 * norm.pdf(grid, 0.2, 0.06) + 0.5 * norm.pdf(grid, 0.55, 0.08); a /= a.sum()
b = norm.pdf(grid, 0.7, 0.1); b /= b.sum()
C = (grid[:, None] - grid[None, :]) ** 2


def sinkhorn(eps, iters=5000):
    K = np.exp(-C / eps)
    v = np.ones(n)
    for _ in range(iters):
        u = a / (K @ v)
        v = b / (K.T @ u)
    return u[:, None] * K * v[None, :]


Fa, Gb = np.cumsum(a), np.cumsum(b)
# exact monotone (north-west corner) coupling for 1-D
P_exact = np.zeros((n, n)); i = j = 0; ra, rb = a.copy(), b.copy()
while i < n and j < n:
    m = min(ra[i], rb[j]); P_exact[i, j] += m; ra[i] -= m; rb[j] -= m
    if ra[i] <= 1e-15: i += 1
    if j < n and rb[j] <= 1e-15: j += 1
W2_exact = np.sqrt((P_exact * C).sum())
print("exact W2", W2_exact)

fig, axs = new_figure(1, 4, figsize=(11.5, 3.4))
for ax, (eps, P) in zip(axs, [(0.05, sinkhorn(0.05)), (0.005, sinkhorn(0.005)), (0.0005, sinkhorn(0.0005)),
                              (None, P_exact)]):
    ax.imshow(P, origin="lower", extent=[0, 1, 0, 1], cmap="Blues", aspect="equal",
              vmax=P.max())
    ax.grid(False)
    cost = np.sqrt((P * C).sum())
    title = (r"$\varepsilon=%g$" % eps) if eps else "exact (monotone)"
    ax.set_title(title + "\n" + r"$\sqrt{\langle \mathbf{C},\mathbf{P}\rangle}=%.4f$" % cost, fontsize=10)
    ax.set_xlabel("target position $y$")
    print(eps, cost)
axs[0].set_ylabel("source position $x$")
save(fig, "ch55-sinkhorn")
