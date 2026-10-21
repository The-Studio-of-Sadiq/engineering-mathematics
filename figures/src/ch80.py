import sys, pathlib; sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
from style import plt, np, PALETTE, INK, MUTED, new_figure, diagram, node, box, arrow, save
from scipy.stats import norm, binom
from matplotlib.patches import Circle, Rectangle


def h2(p):
    p = np.clip(p, 1e-15, 1 - 1e-15)
    return -p * np.log2(p) - (1 - p) * np.log2(1 - p)


# Gauss-Hermite nodes for E[f(Z)], Z ~ N(0,1)
gh_x, gh_w = np.polynomial.hermite_e.hermegauss(120)
gh_w = gh_w / gh_w.sum()


def c_biawgn(snr):
    """Capacity of binary-input (+-1) AWGN with noise variance 1/snr, soft decisions."""
    s2 = 1.0 / snr
    y = 1.0 + np.sqrt(s2) * gh_x
    return 1.0 - np.sum(gh_w * np.logaddexp(0.0, -2.0 * y / s2)) / np.log(2)


# --- Figure 1: capacity curves and the Shannon limit --------------------------------
snr_db = np.linspace(-10, 20, 400)
snr = 10 ** (snr_db / 10)
c_gauss = 0.5 * np.log2(1 + snr)
c_soft = np.array([c_biawgn(s) for s in snr])
c_hard = 1 - h2(norm.sf(np.sqrt(snr)))

fig, (a1, a2) = new_figure(1, 2, figsize=(10.2, 4.3))
a1.plot(snr_db, c_gauss, color=PALETTE[0], label=r"any input: $\frac{1}{2}\log_2(1+\mathrm{SNR})$")
a1.plot(snr_db, c_soft, color=PALETTE[2], label="binary input, soft decisions")
a1.plot(snr_db, c_hard, color=PALETTE[1], ls="--", label="binary input, hard decisions (BSC)")
a1.set_xlabel("SNR per real sample (dB)")
a1.set_ylabel("capacity (bits per channel use)")
a1.set_xlim(-10, 20); a1.set_ylim(0, 3.5)
a1.set_title("Capacity of a Gaussian-noise channel")
a1.legend(loc="upper left", fontsize=8.6)

# right: rate vs Eb/N0
R = np.linspace(0.005, 3.0, 400)
eb = (2 ** (2 * R) - 1) / (2 * R)
a2.fill_betweenx(R, -3, 10 * np.log10(eb), color=PALETTE[1], alpha=0.08)
a2.plot(10 * np.log10(eb), R, color=PALETTE[0], label="Shannon bound (any input)")
ebs = snr / (2 * c_soft)
a2.plot(10 * np.log10(ebs), c_soft, color=PALETTE[2], label="binary-input bound")
a2.axvline(10 * np.log10(np.log(2)), color=INK, ls=":", lw=1.5)
a2.text(-1.45, 2.55, "Shannon limit\n$-1.59$ dB", fontsize=9, color=INK)
a2.text(-2.8, 1.1, "no code\ncan work\nhere", fontsize=9.5, color=PALETTE[1])
a2.plot([9.59], [1.0], "s", color=PALETTE[4], ms=7)
a2.annotate("uncoded BPSK,\nbit error rate $10^{-5}$", xy=(9.59, 1.0), xytext=(6.2, 1.7),
            arrowprops=dict(arrowstyle="->", color=INK), fontsize=9)
a2.annotate("", xy=(1.76, 1.0), xytext=(9.4, 1.0),
            arrowprops=dict(arrowstyle="<->", color=PALETTE[4], lw=1.4))
a2.text(5.6, 0.82, "7.8 dB of room for coding", fontsize=8.8, color=PALETTE[4], ha="center")
a2.set_xlim(-3, 12); a2.set_ylim(0, 3.0)
a2.set_xlabel(r"energy per information bit $E_b/N_0$ (dB)")
a2.set_ylabel("rate $R$ (bits per real channel use)")
a2.set_title("The price of a bit")
a2.legend(loc="lower right", fontsize=8.6)
save(fig, "ch80-capacity-curves")

# --- Figure 2: sphere packing and the concentration of noise ---------------------------
fig, (b1, b2) = new_figure(1, 2, figsize=(10.2, 4.4), gridspec_kw=dict(width_ratios=[1.05, 1]))
b1.set_xlim(0, 10); b1.set_ylim(0, 8); b1.set_aspect("equal"); b1.axis("off")
b1.add_patch(Rectangle((0.2, 0.2), 9.6, 7.0, fc="#f8fafc", ec=INK, lw=1.5))
rng = np.random.default_rng(7)
centres = []
r = 0.78
while len(centres) < 60:
    c = rng.uniform([0.2 + r, 0.2 + r], [9.8 - r, 7.2 - r])
    if all(np.hypot(*(c - d)) > 2 * r + 0.06 for d in centres):
        centres.append(c)
    if len(centres) >= 17:
        break
for i, (x, y) in enumerate(centres):
    hl = i == 0
    b1.add_patch(Circle((x, y), r, fc=PALETTE[0], alpha=0.28 if hl else 0.12, ec=PALETTE[0], lw=1.6 if hl else 1))
    b1.plot([x], [y], "o", color=PALETTE[0], ms=4.5)
x0, y0 = centres[0]
yx, yy = x0 + 0.45, y0 + 0.4
b1.plot([yx], [yy], "*", color=PALETTE[1], ms=11, zorder=5)
b1.annotate("received word $\\mathbf{y}$", xy=(yx, yy), xytext=(yx - 1.0, 7.5), fontsize=9, ha="center",
            arrowprops=dict(arrowstyle="->", color=PALETTE[1]), color=PALETTE[1])
b1.text(0.25, 7.5, r"all $2^n$ words of length $n$", ha="left", fontsize=9.5, color=INK)
b1.text(5.0, -0.25, r"each ball: $\approx 2^{nH_b(\varepsilon)}$ likely noisy versions of one codeword" "\n"
        r"so at most $2^{n(1-H_b(\varepsilon))}$ codewords fit", ha="center", va="top", fontsize=9.2, color=INK)
b1.set_title("Codewords as centres of noise balls")

eps = 0.11
for n, col in zip((20, 200, 2000), (PALETTE[4], PALETTE[2], PALETTE[3])):
    k = np.arange(n + 1)
    pmf = binom.pmf(k, n, eps)
    b2.plot(k / n, pmf * n, color=col, lw=1.8, label=f"$n={n}$")
b2.axvline(eps, color=INK, ls=":", lw=1.3)
b2.text(eps + 0.008, 50, r"$\varepsilon=0.11$", fontsize=9)
b2.set_xlim(0, 0.35); b2.set_ylim(0, 62)
b2.set_xlabel("fraction of bits flipped in a block")
b2.set_ylabel("probability density")
b2.set_title("Long blocks make noise predictable")
b2.legend(loc="upper right")
save(fig, "ch80-sphere-packing")
