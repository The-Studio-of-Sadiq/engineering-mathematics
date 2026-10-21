import sys, pathlib; sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
from style import plt, np, PALETTE, INK, MUTED, new_figure, diagram, node, box, arrow, save
from scipy import stats

th = np.linspace(0, 1, 501)

# --- Figure 1: prior, likelihood, posterior for a beta-binomial coin ---------------
a, b, k, n = 2, 2, 7, 10
prior = stats.beta(a, b).pdf(th)
lik = th ** k * (1 - th) ** (n - k)
lik = lik / np.trapezoid(lik, th)                      # scaled to unit area for display
post = stats.beta(a + k, b + n - k).pdf(th)
fig, ax = new_figure(figsize=(7.0, 3.9))
ax.plot(th, prior, color=PALETTE[6], ls="--", label=r"prior  Beta(2, 2)")
ax.plot(th, lik, color=PALETTE[4], ls=":", lw=2.4, label=r"likelihood  $\theta^7(1-\theta)^3$ (scaled)")
ax.plot(th, post, color=PALETTE[0], lw=2.6, label=r"posterior  Beta(9, 5)")
lo, hi = stats.beta(9, 5).ppf([0.025, 0.975])
m = (th >= lo) & (th <= hi)
ax.fill_between(th[m], post[m], color=PALETTE[0], alpha=0.13, lw=0)
for x, lab, c in [(0.5, "prior mean 0.5", PALETTE[6]), (9 / 14, "posterior mean 0.643", PALETTE[0]),
                  (0.7, "MLE 0.7", PALETTE[4])]:
    ax.axvline(x, color=c, lw=1, alpha=0.6)
ax.text(0.488, 3.45, "prior\nmean 0.5", ha="right", fontsize=8.5, color=PALETTE[6])
ax.text(9 / 14 - 0.012, 3.45, "post.\nmean 0.643", ha="right", fontsize=8.5, color=PALETTE[0])
ax.text(0.71, 3.45, "MLE\n0.7", ha="left", fontsize=8.5, color=PALETTE[4])
ax.text(0.62, 0.6, "95 % credible\n[0.386, 0.861]", ha="center", fontsize=8.5, color=PALETTE[0])
ax.set_xlabel(r"heads probability $\theta$")
ax.set_ylabel("density")
ax.set_ylim(0, 4.0)
ax.set_title("7 heads in 10 tosses: the posterior is a compromise")
ax.legend(loc="upper left", fontsize=8.5)
save(fig, "ch35-prior-posterior")

# --- Figure 2: sequential updating --------------------------------------------------
rng = np.random.default_rng(35)
flips = rng.random(256) < 0.7
fig, ax = new_figure(figsize=(7.0, 3.8))
checkpoints = [0, 1, 4, 16, 64, 256]
cmap = plt.get_cmap("Blues")
for i, nn in enumerate(checkpoints):
    kk = int(flips[:nn].sum())
    d = stats.beta(1 + kk, 1 + nn - kk).pdf(th)
    ax.plot(th, d, color=cmap(0.3 + 0.7 * i / (len(checkpoints) - 1)), lw=2,
            label=f"n = {nn:3d}, heads = {kk}")
ax.axvline(0.7, color=PALETTE[1], ls="--", lw=1.2)
ax.text(0.73, 9.0, r"true $\theta=0.7$", color=PALETTE[1], fontsize=9)
ax.set_xlabel(r"$\theta$")
ax.set_ylabel("posterior density")
ax.set_ylim(0, 14.5)
ax.set_title("Uniform prior, updated flip by flip")
ax.legend(loc="upper left", fontsize=8.5)
save(fig, "ch35-sequential")

# --- Figure 3: Lindley's paradox -----------------------------------------------------
ns = np.unique(np.round(np.logspace(2, 6, 40)).astype(int))
bf, pv = [], []
for nn in ns:
    kk = int(round(nn / 2 + 0.98 * np.sqrt(nn)))
    bf.append(np.exp(stats.binom.logpmf(kk, nn, 0.5) + np.log(nn + 1)))
    pv.append(2 * stats.binom.sf(kk - 1, nn, 0.5))
fig, (a1, a2) = new_figure(1, 2, figsize=(8.4, 3.5))
a1.semilogx(ns, pv, color=PALETTE[4])
a1.axhline(0.05, color=MUTED, ls="--", lw=1)
a1.set_ylim(0, 0.1)
a1.set_xlabel("number of tosses $n$")
a1.set_ylabel("two-sided p-value")
a1.set_title("p-value: always ≈ 0.05")
a2.loglog(ns, bf, color=PALETTE[0])
a2.axhline(1, color=MUTED, ls="--", lw=1)
a2.set_ylim(0.7, 200)
a2.text(2e4, 0.8, "evidence neutral", fontsize=8.5, color=MUTED)
a2.set_xlabel("number of tosses $n$")
a2.set_ylabel(r"Bayes factor $B_{01}$ (fair : not fair)")
a2.set_title("Bayes factor: ever more 'fair'")
save(fig, "ch35-lindley")
