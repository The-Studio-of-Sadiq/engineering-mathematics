"""Figures for Chapter 91 — Inverse Problems: Reconstructing the Invisible."""
import sys, pathlib; sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
from style import plt, np, PALETTE, INK, MUTED, new_figure, diagram, node, box, arrow, save

BLUE, RED, GREEN, PURPLE, ORANGE, TEAL, GREY = PALETTE

# ---------------------------------------------------------------- the test problem (same as the chapter's numbers)
n = 200
x = (np.arange(n) + 0.5) / n
h = 1 / n
s = 0.01
K = h * np.exp(-(x[:, None] - x[None, :]) ** 2 / (2 * s * s)) / (np.sqrt(2 * np.pi) * s)
f = (np.where((x > 0.15) & (x < 0.35), 1.0, 0) + 0.8 * np.maximum(0, 1 - np.abs(x - 0.6) / 0.1)
     + 0.6 * np.exp(-((x - 0.85) / 0.04) ** 2))
rng = np.random.default_rng(0)
g0 = K @ f
noise = 0.01 * np.linalg.norm(g0) / np.sqrt(n) * rng.standard_normal(n)
g = g0 + noise
delta = np.linalg.norm(noise)
U, sv, Vt = np.linalg.svd(K)
beta = U.T @ g


def tik(lam):
    return Vt.T @ (sv / (sv**2 + lam) * beta)


f_naive = np.linalg.solve(K, g)

# ---------------------------------------------------------------- figure 1: deblurring
fig, axs = new_figure(1, 3, figsize=(12.0, 3.9))
ax = axs[0]
ax.plot(x, f, color=INK, lw=2.2, label="true signal $f$")
ax.plot(x, g, color=ORANGE, lw=1.3, label="blurred + 1% noise $g$")
ax.set_ylim(-0.25, 1.35); ax.set_xlabel("$x$")
ax.legend(fontsize=8.8, loc="upper right")
ax.set_title("What the instrument records", fontsize=11)

ax = axs[1]
ax.plot(x, f_naive / 1e5, color=RED, lw=1.0)
ax.set_xlabel("$x$"); ax.set_ylabel(r"naive $\mathbf{K}^{-1}\mathbf{g}$  ($\times10^{5}$)")
ax.set_title("Exact inverse: noise amplified $10^5$-fold", fontsize=11)

ax = axs[2]
ax.plot(x, f, color=INK, lw=2.6, alpha=0.35, label="true $f$")
for lam, col, lab, lw in [(1e-6, PURPLE, r"$\lambda=10^{-6}$ (too small)", 0.8),
                      (5e-3, GREEN, r"$\lambda=5\times10^{-3}$ (about right)", 1.8),
                      (0.3, BLUE, r"$\lambda=0.3$ (too large)", 1.8)]:
    ax.plot(x, tik(lam), color=col, lw=lw, label=lab, alpha=0.55 if lw < 1 else 1)
ax.set_ylim(-0.6, 1.6); ax.set_xlabel("$x$")
ax.legend(fontsize=8.6, loc="upper center", ncol=2, bbox_to_anchor=(0.5, -0.17))
ax.set_title("Tikhonov at three strengths", fontsize=11)
save(fig, "ch91-deblurring")

# ---------------------------------------------------------------- figure 2: L-curve and parameter choice
lams = np.logspace(-10, 1, 221)
res = np.array([np.linalg.norm(K @ tik(l) - g) for l in lams])
nrm = np.array([np.linalg.norm(tik(l)) for l in lams])
err = np.array([np.linalg.norm(tik(l) - f) / np.linalg.norm(f) for l in lams])
t = np.log(lams)
rho, eta = np.log(res), np.log(nrm)
d1r, d1e = np.gradient(rho, t), np.gradient(eta, t)
d2r, d2e = np.gradient(d1r, t), np.gradient(d1e, t)
kap = (d1r * d2e - d2r * d1e) / (d1r**2 + d1e**2) ** 1.5
i_corner = np.argmax(kap)
i_disc = np.argmin(np.abs(res - delta))
i_best = np.argmin(err)

fig, axs = new_figure(1, 2, figsize=(10.6, 4.2))
ax = axs[0]
ax.loglog(res, nrm, color=BLUE, lw=2)
for i, col, lab, off in [(i_corner, PURPLE, "L-curve corner", (8, 32)),
                         (i_disc, GREEN, "discrepancy principle", (12, -30)),
                         (i_best, RED, "best possible", (12, -48))]:
    ax.plot(res[i], nrm[i], "o", color=col, ms=7, zorder=4)
    ax.annotate(lab + rf"  ($\lambda\approx${lams[i]:.0e})", xy=(res[i], nrm[i]), xytext=off,
                textcoords="offset points", fontsize=9, color=col)
for lam in [1e-8, 1e-1, 1]:
    j = np.argmin(np.abs(lams - lam))
    ax.plot(res[j], nrm[j], ".", color=MUTED, ms=6)
    ax.annotate(rf"$10^{{{int(np.log10(lam))}}}$", xy=(res[j], nrm[j]), xytext=(5, 3),
                textcoords="offset points", fontsize=8.5, color=MUTED)
ax.axvline(delta, color=GREEN, lw=1, ls=":")
ax.text(delta * 1.05, 40, r"noise level $\delta$", color=GREEN, fontsize=9, rotation=90, va="bottom")
ax.set_xlabel(r"residual $\Vert\mathbf{K}\mathbf{f}_\lambda-\mathbf{g}\Vert$")
ax.set_ylabel(r"solution size $\Vert\mathbf{f}_\lambda\Vert$")
ax.set_title("The L-curve", fontsize=11)

ax = axs[1]
ax.loglog(lams, err, color=INK, lw=2)
for i, col in [(i_corner, PURPLE), (i_disc, GREEN), (i_best, RED)]:
    ax.plot(lams[i], err[i], "o", color=col, ms=7, zorder=4)
ax.set_ylim(0.07, 20)
ax.set_xlabel(r"regularization parameter $\lambda$"); ax.set_ylabel("relative error")
ax.text(1e-9, 3.0, "noise\ndominates", fontsize=9.5, color=MUTED)
ax.text(1e-1, 0.5, "bias\ndominates", fontsize=9.5, color=MUTED)
ax.set_title("Error against $\\lambda$ (needs the true $f$)", fontsize=11)
save(fig, "ch91-l-curve")

# ---------------------------------------------------------------- figure 3: Picard plot
fig, ax = new_figure(figsize=(7.4, 4.2))
k = np.arange(1, n + 1)
ax.semilogy(k, sv, color=BLUE, lw=2, label=r"singular values $\sigma_k$")
ax.semilogy(k, np.abs(beta), ".", color=ORANGE, ms=4, label=r"data coefficients $|\mathbf{u}_k^\mathsf{T}\mathbf{g}|$")
ax.semilogy(k, np.abs(beta) / sv, ".", color=RED, ms=4, label=r"naive solution coefficients $|\mathbf{u}_k^\mathsf{T}\mathbf{g}|/\sigma_k$")
ax.axhline(delta / np.sqrt(n), color=GREEN, lw=1.2, ls="--")
ax.text(178, delta / np.sqrt(n) * 30, "noise floor", color=GREEN, fontsize=9.5, ha="center",
        bbox=dict(fc="white", ec="none", pad=1))
ax.set_xlabel("index $k$"); ax.set_ylim(1e-10, 1e7)
ax.legend(fontsize=9, loc="upper left")
ax.set_title("Picard plot: where the data stop carrying information", fontsize=11)
save(fig, "ch91-picard")
