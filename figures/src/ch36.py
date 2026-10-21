import sys, pathlib; sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
from style import plt, np, PALETTE, INK, MUTED, new_figure, diagram, node, box, arrow, save

# --- Figure 1: Poisson arrivals and the counting process -----------------------------
rng = np.random.default_rng(46)
lam, T = 0.1, 120.0                       # buses per minute, horizon in minutes
gaps = rng.exponential(1 / lam, 40)
arr = np.cumsum(gaps)
arr = arr[arr < T]
fig, (a1, a2) = new_figure(2, 1, figsize=(7.4, 4.6), sharex=True,
                           gridspec_kw=dict(height_ratios=[1, 2.6]))
a1.eventplot(arr, lineoffsets=0, linelengths=0.8, colors=PALETTE[0], lw=2)
a1.set_yticks([])
a1.set_ylim(-0.7, 0.9)
a1.set_title(r"Bus arrivals: a Poisson process, $\lambda=0.1$ per minute")
# highlight the gap that contains a passenger arriving at t = 60
t0 = 60.0
i = np.searchsorted(arr, t0)
g0, g1 = (arr[i - 1] if i > 0 else 0.0), arr[i]
a1.axvspan(g0, g1, color=PALETTE[4], alpha=0.18, lw=0)
a1.plot([t0], [0], marker="v", color=PALETTE[1], ms=9)
a1.text(t0, 0.45, "you arrive", ha="center", fontsize=8.5, color=PALETTE[1])
a1.grid(False)
for s in ["left"]:
    a1.spines[s].set_visible(False)
tt = np.concatenate([[0], arr, [T]])
NN = np.concatenate([[0], np.arange(1, len(arr) + 1), [len(arr)]])
a2.step(tt, NN, where="post", color=PALETTE[0], label=r"$N(t)$")
a2.plot([0, T], [0, lam * T], color=MUTED, ls="--", lw=1.4, label=r"mean $\lambda t$")
a2.axvspan(g0, g1, color=PALETTE[4], alpha=0.18, lw=0)
a2.text((g0 + g1) / 2, 8.6, f"gap holding you:\n{g1 - g0:.0f} min", ha="center", fontsize=8.5, color=PALETTE[4])
a2.set_xlabel("time $t$ (minutes)")
a2.set_ylabel("buses so far")
a2.set_xlim(0, T)
a2.legend(loc="upper left")
save(fig, "ch36-poisson")

# --- Figure 2: Brownian sample paths with the sqrt(t) envelope ----------------------
rng = np.random.default_rng(7)
n, Tb = 2000, 1.0
t = np.linspace(0, Tb, n + 1)
dt = Tb / n
fig, ax = new_figure(figsize=(7.0, 4.0))
for k in range(8):
    W = np.concatenate([[0], np.cumsum(rng.normal(0, np.sqrt(dt), n))])
    ax.plot(t, W, lw=1.0, color=PALETTE[k % 6], alpha=0.85)
for c, ls, lab in [(1, "--", r"$\pm\sqrt{t}$ (68 %)"), (2, ":", r"$\pm2\sqrt{t}$ (95 %)")]:
    ax.plot(t, c * np.sqrt(t), color=INK, ls=ls, lw=1.8, label=lab)
    ax.plot(t, -c * np.sqrt(t), color=INK, ls=ls, lw=1.8)
ax.set_xlabel("time $t$")
ax.set_ylabel("$W(t)$")
ax.set_xlim(0, 1)
ax.set_title("Eight Brownian paths: spread grows like the square root of time")
ax.legend(loc="lower left")
save(fig, "ch36-brownian")

# --- Figure 3: quadratic variation and the Ito correction ---------------------------
rng = np.random.default_rng(3)
N = 2 ** 16
dW = rng.normal(0, np.sqrt(1 / N), N)
W = np.concatenate([[0], np.cumsum(dW)])
levels = np.arange(2, 17)
qv, left, mid = [], [], []
for L in levels:
    step = N // 2 ** L
    Wk = W[::step]
    inc = np.diff(Wk)
    qv.append(np.sum(inc ** 2))
    left.append(np.sum(Wk[:-1] * inc))
    mid.append(np.sum(0.5 * (Wk[:-1] + Wk[1:]) * inc))
mesh = 2.0 ** levels
fig, (a1, a2) = new_figure(1, 2, figsize=(8.6, 3.6))
a1.semilogx(mesh, qv, "o-", ms=3.5, color=PALETTE[0], label=r"$\sum(\Delta W)^2$")
tv = [np.sum(np.abs(np.diff(W[::N // 2 ** L]))) for L in levels]
a1.semilogx(mesh, tv, "s-", ms=3.5, color=PALETTE[1], label=r"$\sum|\Delta W|$")
a1.axhline(1, color=MUTED, ls="--", lw=1.2)
a1.text(3e3, 0.7, "t = 1", color=MUTED, fontsize=9)
a1.set_yscale("log")
a1.set_xlabel("number of subintervals of [0, 1]")
a1.set_title("Squares settle, lengths explode")
a1.legend(loc="upper left", fontsize=9)
WT = W[-1]
a2.semilogx(mesh, left, "o-", ms=3.5, color=PALETTE[0], label="left-point sums (Itô)")
a2.semilogx(mesh, mid, "s-", ms=3.5, color=PALETTE[4], label="trapezoid sums (Stratonovich)")
a2.axhline((WT ** 2 - 1) / 2, color=PALETTE[0], ls="--", lw=1.2)
a2.axhline(WT ** 2 / 2, color=PALETTE[4], ls="--", lw=1.2)
a2.text(mesh[-1], (WT ** 2 - 1) / 2 + 0.04, r"$(W^2-t)/2$", ha="right", va="bottom", color=PALETTE[0], fontsize=9)
a2.text(mesh[-1], WT ** 2 / 2 - 0.02, r"$W^2/2$", ha="right", va="top", color=PALETTE[4], fontsize=9)
a2.set_xlabel("number of subintervals of [0, 1]")
a2.set_title(r"Two answers for $\int_0^1 W\,dW$")
a2.legend(loc="center right", fontsize=8.5)
a2.set_ylim(0.35, 0.98)
save(fig, "ch36-ito")
print("W(1) =", WT)
