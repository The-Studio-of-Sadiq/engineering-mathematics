import sys, pathlib; sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
from style import plt, np, PALETTE, INK, MUTED, new_figure, diagram, node, box, arrow, save
from matplotlib.colors import ListedColormap

rng = np.random.default_rng(88)


def stabilise(z):
    """Topple every site with >= 4 grains until stable; grains falling off the edge are lost.
    Returns the stable grid and the total number of topplings."""
    total = 0
    while True:
        t = z // 4
        if not t.any():
            return z, total
        total += int(t.sum())
        z = z - 4 * t
        z[1:, :] += t[:-1, :]
        z[:-1, :] += t[1:, :]
        z[:, 1:] += t[:, :-1]
        z[:, :-1] += t[:, 1:]


# --- Figure 1: sandpile from a single pile + avalanche-size distribution ---------------
n = 181
z = np.zeros((n, n), dtype=np.int64)
z[n // 2, n // 2] = 2 ** 15
z, _ = stabilise(z)
rows = np.nonzero(z.any(1))[0]
c0, c1 = rows.min() - 2, rows.max() + 3
pattern = z[c0:c1, c0:c1]

# slowly driven pile: drop grains one at a time at random sites
L = 50
g = rng.integers(0, 4, (L, L)).astype(np.int64)
sizes = []
for step in range(60000):
    i, j = rng.integers(0, L, 2)
    g[i, j] += 1
    if g[i, j] >= 4:
        g, s = stabilise(g)
        if step > 10000:
            sizes.append(s)
sizes = np.array(sizes)

fig, (a1, a2) = new_figure(1, 2, figsize=(10.5, 4.6))
cmap = ListedColormap(["#f8fafc", "#bfdbfe", PALETTE[0], "#1e3a8a"])
im = a1.imshow(pattern, cmap=cmap, interpolation="nearest", vmin=-0.5, vmax=3.5)
a1.set_title(r"$2^{15}$ grains dropped on one site, then relaxed")
a1.grid(False); a1.set_xticks([]); a1.set_yticks([])
cb = fig.colorbar(im, ax=a1, ticks=[0, 1, 2, 3], shrink=0.8)
cb.set_label("grains on the site")

bins = np.unique(np.round(np.logspace(0, np.log10(sizes.max() + 1), 28)))
h, e = np.histogram(sizes, bins=bins)
cen = np.sqrt(e[:-1] * e[1:])
ok = h > 0
dens = h[ok] / np.diff(e)[ok] / len(sizes)
a2.loglog(cen[ok], dens, "o", color=PALETTE[1], ms=5, label="simulated avalanches")
fit = (cen[ok] > 3) & (cen[ok] < 300)
slope, icpt = np.polyfit(np.log(cen[ok][fit]), np.log(dens[fit]), 1)
xx = np.logspace(0.3, 3, 20)
a2.loglog(xx, np.exp(icpt) * xx ** slope, "--", color=INK, lw=1.3, label=f"slope {slope:.2f}")
a2.set_xlabel("avalanche size $s$ (number of topplings)")
a2.set_ylabel("probability density $P(s)$")
a2.set_title(f"Avalanches on a {L}×{L} pile: no typical size")
a2.legend(loc="lower left")
save(fig, "ch88-sandpile")
print("avalanches", len(sizes), "max", sizes.max(), "slope", slope, "mean density", g.mean(),
      "pattern radius", (c1 - c0) / 2)

# --- Figure 2: feedback loops --------------------------------------------------------
fig, ax = diagram(figsize=(10.0, 4.6))
ax.set_xlim(-0.5, 14.9)
ax.set_ylim(0, 7.0)


def loop(cx, title, labels, signs, color, sym, word):
    ang = np.deg2rad([90, 330, 210])
    R = 2.3
    pts = [(cx + 1.25 * R * np.cos(a), 3.4 + R * np.sin(a)) for a in ang]
    for (x, y), lab in zip(pts, labels):
        box(ax, x, y, lab, w=2.3, h=0.95, color=color, fontsize=9.3)
    (tx, ty), (rx, ry), (lx, ly) = pts
    hw, hh = 1.15, 0.475
    ends = [((tx + 0.6, ty - hh), (rx, ry + hh)),            # top -> lower right
            ((rx - hw, ry - 0.15), (lx + hw, ly - 0.15)),    # lower right -> lower left
            ((lx, ly + hh), (tx - 0.6, ty - hh))]            # lower left -> top
    for k in range(3):
        (x0, y0), (x1, y1) = pts[k], pts[(k + 1) % 3]
        ax.annotate("", xy=ends[k][1], xytext=ends[k][0], zorder=2,
                    arrowprops=dict(arrowstyle="-|>,head_length=0.6,head_width=0.3", color=color, lw=1.8,
                                    shrinkA=2, shrinkB=2, connectionstyle="arc3,rad=-0.25"))
        mx, my = (x0 + x1) / 2, (y0 + y1) / 2
        dx, dy = mx - cx, my - 3.4
        d = np.hypot(dx, dy)
        ax.text(mx + 0.75 * dx / d, my + 0.75 * dy / d, signs[k], ha="center", va="center",
                fontsize=14, fontweight="bold", color=color)
    ax.add_patch(plt.Circle((cx, 3.4), 0.62, fc="white", ec=color, lw=2))
    ax.text(cx, 3.55, sym, ha="center", va="center", fontsize=16, fontweight="bold", color=color)
    ax.text(cx, 3.12, word, ha="center", va="center", fontsize=8.5, color=color)
    ax.text(cx, 6.75, title, ha="center", va="center", fontsize=11, fontweight="bold", color=INK)


loop(3.3, "Reinforcing loop: a bank run",
     ["Withdrawals", "Bank's\nreserves", "Depositors'\nfear"], ["−", "−", "+"],
     PALETTE[1], "R", "amplifies")
loop(11.1, "Balancing loop: a sandpile",
     ["Slope of\nthe pile", "Avalanche\nsize", "Grains lost\nat the edge"], ["+", "+", "−"],
     PALETTE[0], "B", "stabilises")
ax.text(7.2, 0.1, "Loop polarity = product of the signs on its arrows:  (−)(−)(+) = +  reinforces;"
        "  (+)(+)(−) = −  balances", ha="center", fontsize=9.5, color=MUTED)
save(fig, "ch88-feedback-loops")


# --- Figure 3: Vicsek flocking at low and high noise ---------------------------------
def vicsek(eta, N=400, Lb=10.0, v=0.03, r=1.0, steps=600):
    x = rng.uniform(0, Lb, (N, 2))
    th = rng.uniform(-np.pi, np.pi, N)
    for _ in range(steps):
        d = x[:, None, :] - x[None, :, :]
        d -= Lb * np.round(d / Lb)
        near = (d ** 2).sum(-1) < r * r
        s = near @ np.sin(th); c = near @ np.cos(th)
        th = np.arctan2(s, c) + eta * rng.uniform(-0.5, 0.5, N)
        x = (x + v * np.c_[np.cos(th), np.sin(th)]) % Lb
    va = np.hypot(np.cos(th).mean(), np.sin(th).mean())
    return x, th, va


fig, axes = new_figure(1, 2, figsize=(10.0, 5.0))
for ax, eta, col in zip(axes, [0.3, 5.0], [PALETTE[2], PALETTE[4]]):
    x, th, va = vicsek(eta)
    ax.quiver(x[:, 0], x[:, 1], np.cos(th), np.sin(th), color=col, scale=28, width=0.004,
              headwidth=4, headlength=5)
    ax.set_xlim(0, 10); ax.set_ylim(0, 10); ax.set_aspect("equal")
    ax.set_xticks([]); ax.set_yticks([]); ax.grid(False)
    for s_ in ax.spines.values():
        s_.set_visible(True); s_.set_color(MUTED)
    ax.set_title(f"noise η = {eta}: alignment $v_a$ = {va:.2f}")
    print("vicsek", eta, va)
fig.suptitle("400 self-propelled particles, each copying its neighbours' average heading (plus noise)",
             fontsize=10.5, color=INK)
save(fig, "ch88-vicsek")
