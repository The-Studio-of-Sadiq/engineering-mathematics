"""Chapter 38 figures: estimating pi with darts, Metropolis chains on a 2-D density."""
import sys, pathlib; sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
from style import plt, np, PALETTE, INK, MUTED, new_figure, diagram, node, box, arrow, save

BLUE, RED, GREEN, PURPLE, ORANGE = PALETTE[:5]

# ---------------------------------------------------------------- 1. pi from darts
rng = np.random.default_rng(38)
fig, (ax1, ax2) = new_figure(1, 2, figsize=(10.0, 4.3), gridspec_kw={"width_ratios": [1, 1.35]})
n_show = 1500
p = rng.random((n_show, 2)) * 2 - 1
inside = (p ** 2).sum(1) <= 1
ax1.scatter(*p[inside].T, s=5, color=BLUE, label=f"inside ({inside.sum()})")
ax1.scatter(*p[~inside].T, s=5, color=RED, label=f"outside ({(~inside).sum()})")
t = np.linspace(0, 2 * np.pi, 400)
ax1.plot(np.cos(t), np.sin(t), color=INK, lw=1.5)
ax1.plot([-1, 1, 1, -1, -1], [-1, -1, 1, 1, -1], color=INK, lw=1.5)
ax1.set_aspect("equal"); ax1.set_xlim(-1.05, 1.05); ax1.set_ylim(-1.05, 1.32)
ax1.set_xticks([-1, 0, 1]); ax1.set_yticks([-1, 0, 1])
ax1.set_title(f"{n_show} darts:  4 x {inside.sum()}/{n_show} = {4 * inside.mean():.3f}")
ax1.legend(loc="upper center", ncol=2, fontsize=9, markerscale=3, handletextpad=0.2)

N = 100_000
hits = ((rng.random((N, 2)) * 2 - 1) ** 2).sum(1) <= 1
n = np.arange(1, N + 1)
est = 4 * np.cumsum(hits) / n
sd = 4 * np.sqrt(np.pi / 4 * (1 - np.pi / 4))
ax2.fill_between(n, np.pi - 2 * sd / np.sqrt(n), np.pi + 2 * sd / np.sqrt(n), color="#e5edff", label=r"$\pi \pm 2\sigma/\sqrt{N}$")
ax2.plot(n, est, color=BLUE, lw=1.4, label="running estimate")
ax2.axhline(np.pi, color=RED, lw=1.2, ls="--", label=r"$\pi$")
ax2.set_xscale("log"); ax2.set_ylim(2.4, 3.9); ax2.set_xlim(10, N)
ax2.set_xlabel("number of darts N (log scale)"); ax2.set_ylabel("estimate of pi")
ax2.set_title("The error shrinks like 1/sqrt(N)")
ax2.legend(loc="upper right", fontsize=9)
save(fig, "ch38-pi-darts")

# ---------------------------------------------------------------- 2. Metropolis chains
def logpi(x, y):  # a curved "banana" density
    return -x ** 2 / 8 - (y - 0.25 * x ** 2) ** 2 / (2 * 0.5 ** 2)

gx, gy = np.meshgrid(np.linspace(-6, 6, 300), np.linspace(-2.5, 9, 300))
dens = np.exp(logpi(gx, gy))

fig, axes = new_figure(1, 3, figsize=(11.0, 4.1), sharey=True)
steps, labels = [0.08, 2.5, 25.0], ["step too small", "step about right", "step too large"]
for ax, s, lab, c in zip(axes, steps, labels, [ORANGE, GREEN, PURPLE]):
    r = np.random.default_rng(7)
    x = np.array([4.5, -2.0]); lp = logpi(*x); chain = [x]; acc = 0
    for _ in range(3000):
        prop = x + s * r.normal(size=2)
        lq = logpi(*prop)
        if np.log(r.random()) < lq - lp:
            x, lp = prop, lq; acc += 1
        chain.append(x)
    ch = np.array(chain)
    ax.contourf(gx, gy, dens, levels=8, cmap="Blues", alpha=0.55)
    ax.plot(ch[:, 0], ch[:, 1], color=c, lw=0.6, alpha=0.8)
    ax.scatter(ch[:, 0], ch[:, 1], s=2, color=c, zorder=3)
    ax.plot(*ch[0], marker="*", color=RED, ms=12, zorder=4)
    ax.set_title(f"{lab} (s = {s})\nacceptance {acc / 3000:.1%}", fontsize=11)
    ax.set_xlabel("x"); ax.set_xlim(-6, 6); ax.set_ylim(-2.5, 9)
axes[0].set_ylabel("y")
save(fig, "ch38-metropolis-chains")
