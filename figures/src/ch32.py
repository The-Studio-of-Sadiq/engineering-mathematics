import sys, pathlib; sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
from style import plt, np, PALETTE, INK, MUTED, new_figure, diagram, node, box, arrow, save
from scipy import stats

# --- Figure 1: tail bounds against the true binomial tail -----------------------
n, p = 100, 0.5
mu, var = n * p, n * p * (1 - p)
k = np.arange(52, 96)
a = k / n
true = stats.binom.sf(k - 1, n, p)
markov = np.minimum(1, mu / k)
cheb = np.minimum(1, var / (k - mu) ** 2)
kl = a * np.log(a / p) + (1 - a) * np.log((1 - a) / (1 - p))
chern = np.exp(-n * kl)
fig, ax = new_figure(figsize=(7.4, 4.4))
ax.semilogy(k, markov, color=PALETTE[6], label=r"Markov  $\mu/k$")
ax.semilogy(k, cheb, color=PALETTE[4], label=r"Chebyshev  $\sigma^2/(k-\mu)^2$")
ax.semilogy(k, chern, color=PALETTE[0], label=r"Chernoff  $e^{-n\,D(k/n\,\Vert\,p)}$")
ax.semilogy(k, true, color=PALETTE[1], ls="none", marker="o", ms=3.5, label=r"true  $P(X\geq k)$")
ax.set_ylim(1e-20, 2)
ax.set_xlabel(r"threshold $k$   ($X\sim$ Binomial$(100,\ 1/2)$, mean 50)")
ax.set_ylabel("probability (log scale)")
ax.set_title("Three bounds on the same tail")
ax.annotate("polynomial: useless\nfar out", xy=(85, cheb[85 - 52]), xytext=(78, 1e-6),
            arrowprops=dict(arrowstyle="->", color=MUTED), color=MUTED, fontsize=9)
ax.annotate("exponential: tracks the\ntrue tail's slope", xy=(80, chern[80 - 52]), xytext=(56, 1e-11),
            arrowprops=dict(arrowstyle="->", color=MUTED), color=MUTED, fontsize=9)
ax.legend(loc="lower left")
save(fig, "ch32-tail-bounds")

# --- Figure 2: the birthday curve ------------------------------------------------
M = 365
kk = np.arange(1, 81)
exact = 1 - np.cumprod((M - np.arange(80)) / M)
approx = 1 - np.exp(-kk * (kk - 1) / (2 * M))
fig, ax = new_figure(figsize=(7.4, 4.0))
ax.plot(kk, exact, color=PALETTE[0], label=r"exact  $1-\prod_{i<k}(1-i/365)$")
ax.plot(kk, approx, color=PALETTE[1], ls="--", label=r"approximation  $1-e^{-k(k-1)/730}$")
ax.axhline(0.5, color=MUTED, lw=1, ls=":")
ax.plot([23], [exact[22]], "o", color=INK, zorder=5)
ax.annotate(f"k = 23: {exact[22]:.3f}", xy=(23, exact[22]), xytext=(30, 0.32),
            arrowprops=dict(arrowstyle="->", color=INK), fontsize=9.5)
ax.axvline(np.sqrt(M), color=PALETTE[2], lw=1.2, ls="-.")
ax.text(np.sqrt(M) - 1, 0.86, r"$\sqrt{365}\approx 19$", color=PALETTE[2], ha="right", fontsize=9.5)
ax.set_xlabel("number of people $k$"); ax.set_ylabel("P(some shared birthday)")
ax.set_xlim(0, 80); ax.set_ylim(0, 1.03)
ax.set_title("Collisions arrive at the square root")
ax.legend(loc="lower right")
save(fig, "ch32-birthday-curve")
