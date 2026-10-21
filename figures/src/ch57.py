import sys, pathlib; sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
from style import plt, np, PALETTE, INK, MUTED, new_figure, diagram, node, box, arrow, save
from matplotlib.patches import FancyArrowPatch, Circle

# ---- figure 1: the modular clock, 7 + 8 = 3 (mod 12), with the residue classes on spokes
fig, ax = diagram(figsize=(7.6, 5.2))
ax.set_xlim(-4.6, 6.4); ax.set_ylim(-3.6, 3.9)
R = 2.6
ax.add_patch(Circle((0, 0), R, fc="#f8fafc", ec=INK, lw=2))


def pos(h, r=R):
    t = np.pi / 2 - 2 * np.pi * h / 12
    return r * np.cos(t), r * np.sin(t)


for h in range(12):
    x, y = pos(h, 0.82 * R)
    ax.text(x, y, str(h), ha="center", va="center", fontsize=13, color=INK, weight="bold")
    x0, y0 = pos(h, 0.95 * R); x1, y1 = pos(h, R)
    ax.plot([x0, x1], [y0, y1], color=INK, lw=1.5)
# arc: start at 7, add 8 hours, land on 3
th = np.linspace(np.pi / 2 - 2 * np.pi * 7 / 12, np.pi / 2 - 2 * np.pi * 15 / 12, 200)
ax.plot(1.12 * R * np.cos(th[:-6]), 1.12 * R * np.sin(th[:-6]), color=PALETTE[1], lw=2.4)
ax.annotate("", xy=(1.12 * R * np.cos(th[-1]), 1.12 * R * np.sin(th[-1])),
            xytext=(1.12 * R * np.cos(th[-7]), 1.12 * R * np.sin(th[-7])),
            arrowprops=dict(arrowstyle="-|>", color=PALETTE[1], lw=2.4, mutation_scale=16))
for h, c in [(7, PALETTE[0]), (3, PALETTE[1])]:
    x, y = pos(h, 0.82 * R)
    ax.add_patch(Circle((x, y), 0.33, fc="none", ec=c, lw=2.2))
ax.text(0, 0.25, "7 + 8 ≡ 3", ha="center", fontsize=13, color=PALETTE[1], weight="bold")
ax.text(0, -0.35, "(mod 12)", ha="center", fontsize=11, color=PALETTE[1])
ax.text(-1.9, 3.45, "+8 hours", fontsize=10, color=PALETTE[1])
# residue classes listed on the right
ax.text(3.35, 3.3, "each hour is a residue class:", fontsize=10, color=INK, weight="bold")
rows = [(0, "…, −12, 0, 12, 24, …"), (3, "…, −9, 3, 15, 27, …"), (7, "…, −5, 7, 19, 31, …"),
        (11, "…, −1, 11, 23, 35, …")]
for i, (h, s) in enumerate(rows):
    ax.text(3.35, 2.6 - 0.62 * i, f"[{h}] = {s}", fontsize=10, color=PALETTE[0] if h != 3 else PALETTE[1])
ax.text(3.35, -0.25, "multiplication also wraps:", fontsize=10, color=INK, weight="bold")
ax.text(3.35, -0.85, "5 × 5 = 25 ≡ 1   (5 is its own inverse)", fontsize=10, color=PALETTE[2])
ax.text(3.35, -1.45, "3 × 4 = 12 ≡ 0   (zero divisors!)", fontsize=10, color=PALETTE[1])
ax.text(3.35, -2.05, "units mod 12: {1, 5, 7, 11}", fontsize=10, color=PALETTE[2])
ax.text(3.35, -2.65, "φ(12) = 4", fontsize=10, color=PALETTE[2])
save(fig, "ch57-modular-clock")

# ---- figure 2: multiplication-mod-n "times table" circles
fig, axes = plt.subplots(1, 3, figsize=(10.5, 3.9), constrained_layout=True)
for ax, (k, n) in zip(axes, [(2, 200), (3, 200), (34, 200)]):
    t = 2 * np.pi * np.arange(n) / n
    x, y = np.sin(t), np.cos(t)
    for i in range(n):
        j = (k * i) % n
        ax.plot([x[i], x[j]], [y[i], y[j]], color=PALETTE[0], lw=0.45, alpha=0.75)
    ax.plot(np.sin(np.linspace(0, 2 * np.pi, 300)), np.cos(np.linspace(0, 2 * np.pi, 300)), color=INK, lw=1)
    ax.set_aspect("equal"); ax.axis("off")
    ax.set_title(f"i → {k}·i mod {n}", fontsize=11, color=INK)
fig.suptitle("Join point i to point k·i (mod n) on a circle of n points", fontsize=11.5, color=INK, weight="bold")
save(fig, "ch57-times-table")

# ---- figure 3: prime gaps
N = 200000
sieve = np.ones(N + 1, bool); sieve[:2] = False
for p in range(2, int(N ** 0.5) + 1):
    if sieve[p]:
        sieve[p * p::p] = False
P = np.nonzero(sieve)[0]
gaps = np.diff(P)
fig, ax = new_figure(figsize=(7.6, 4.2))
ax.scatter(P[:-1], gaps, s=4, color=PALETTE[0], alpha=0.35, lw=0, label="gap $p_{k+1}-p_k$")
xs = np.geomspace(3, N, 400)
ax.plot(xs, np.log(xs), color=PALETTE[2], lw=2.2, label=r"average gap $\approx \ln p$")
ax.plot(xs, np.log(xs) ** 2, color=PALETTE[1], lw=2, ls="--", label=r"Cramér scale $(\ln p)^2$")
rec = [(P[i], gaps[i]) for i in range(len(gaps)) if gaps[i] == gaps[:i + 1].max() and (i == 0 or gaps[i] > gaps[:i].max())]
rx, ry = zip(*rec)
ax.plot(rx, ry, "o", color=PALETTE[4], ms=5, label="record gaps")
ax.set_xscale("log"); ax.set_xlim(2, N); ax.set_ylim(0, 160)
ax.set_xlabel("prime $p$ (log scale)"); ax.set_ylabel("distance to next prime")
ax.set_title("Prime gaps below 200,000: small on average, wild in detail")
ax.legend(loc="upper left")
save(fig, "ch57-prime-gaps")
