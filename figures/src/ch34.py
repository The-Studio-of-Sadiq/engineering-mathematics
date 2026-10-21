import sys, pathlib; sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
from style import plt, np, PALETTE, INK, MUTED, new_figure, diagram, node, box, arrow, save

# --- Figure 1: weather chain state diagram -------------------------------------------
P = np.array([[0.7, 0.2, 0.1],
              [0.3, 0.4, 0.3],
              [0.2, 0.3, 0.5]])
names = ["Sunny", "Cloudy", "Rainy"]
pos = {0: (2.2, 3.9), 1: (7.8, 3.9), 2: (5.0, 1.0)}
cols = [PALETTE[4], PALETTE[6], PALETTE[0]]

fig, ax = diagram(figsize=(7.0, 4.4))
ax.set_ylim(-0.6, 5.6)
R = 0.75
for i in range(3):
    x, y = pos[i]
    node(ax, x, y, names[i], color=cols[i], r=R, fontsize=10.5, lw=2.4)

# directed edges between distinct states, curved so the two directions separate
for i in range(3):
    for j in range(3):
        if i == j:
            continue
        x0, y0 = pos[i]
        x1, y1 = pos[j]
        th = np.arctan2(y1 - y0, x1 - x0)
        a0, a1_ = th - 0.32, th + np.pi + 0.32      # leave/enter the circles off-axis
        p0 = (x0 + R * np.cos(a0), y0 + R * np.sin(a0))
        p1 = (x1 + R * np.cos(a1_), y1 + R * np.sin(a1_))
        ax.annotate("", xy=p1, xytext=p0, zorder=4,
                    arrowprops=dict(arrowstyle="-|>,head_length=0.7,head_width=0.35", color=cols[i], lw=1.6,
                                    shrinkA=1, shrinkB=1, connectionstyle="arc3,rad=0.18"))
        # label near the curved arc: offset midpoint perpendicular to the chord
        mx, my = (x0 + x1) / 2, (y0 + y1) / 2
        dx, dy = x1 - x0, y1 - y0
        L = np.hypot(dx, dy)
        nx, ny = dy / L, -dx / L          # right-hand normal (arc3 positive rad bends this way)
        off = 0.22 * L * 0.5 + 0.12
        ax.text(mx + nx * off, my + ny * off, f"{P[i, j]:.1f}", ha="center", va="center",
                fontsize=9.5, color=cols[i], fontweight="bold",
                bbox=dict(fc="white", ec="none", pad=0.6), zorder=5)

# self-loops
loop = {0: (-1, 1), 1: (1, 1), 2: (0, -1)}
for i in range(3):
    x, y = pos[i]
    sx, sy = loop[i]
    if sy > 0:
        cx, cy = x + 0.55 * sx, y + 0.85
    else:
        cx, cy = x, y - 1.05
    ax.add_patch(plt.Circle((cx, cy), 0.38, fc="none", ec=cols[i], lw=1.6, zorder=1))
    tx = cx + (0.65 * sx if sx else 0.75)
    ty = cy + (0.25 if sy > 0 else 0.0)
    ax.text(tx, ty, f"{P[i, i]:.1f}", ha="center", va="center", fontsize=9.5,
            color=cols[i], fontweight="bold")

ax.text(5.0, 5.45, "Tomorrow's weather depends only on today's", ha="center", fontsize=10, color=MUTED)
save(fig, "ch34-weather-chain")

# --- Figure 2: convergence to the stationary distribution ---------------------------
pi = np.array([21, 13, 12]) / 46
n = np.arange(0, 16)
fig, (a1, a2) = new_figure(1, 2, figsize=(8.4, 3.6))
starts = [np.eye(3)[k] for k in range(3)]
for k, mu in enumerate(starts):
    traj, tv = [], []
    m = mu.copy()
    for _ in n:
        traj.append(m[0])
        tv.append(0.5 * np.abs(m - pi).sum())
        m = m @ P
    a1.plot(n, traj, "o-", ms=3.5, color=cols[k], label=f"start {names[k]}")
    a2.semilogy(n, tv, "o-", ms=3.5, color=cols[k], label=f"start {names[k]}")
a1.axhline(pi[0], color=INK, ls="--", lw=1.2)
a1.text(15, pi[0] - 0.07, r"$\pi_{\rm Sunny}=21/46$", ha="right", fontsize=9.5, color=INK)
a1.set_xlabel("day $n$")
a1.set_ylabel("P(Sunny on day $n$)")
a1.set_title("Every start forgets itself")
a1.set_ylim(-0.03, 1.03)
a1.legend(loc="upper right", fontsize=8.5)
lam2 = (3 + np.sqrt(3)) / 10
a2.semilogy(n, 0.74 * lam2 ** n, color=MUTED, ls=":", lw=1.6, label=r"slope $\lambda_2^n$, $\lambda_2\approx0.473$")
a2.set_xlabel("day $n$")
a2.set_ylabel(r"distance to $\pi$ (total variation)")
a2.set_title("…at a geometric rate")
a2.legend(fontsize=8.5, loc="lower left")
save(fig, "ch34-convergence")

# --- Figure 3: gambler's ruin sample paths -----------------------------------------
rng = np.random.default_rng(34)
Ncap, k0 = 10, 4
fig, ax = new_figure(figsize=(7.0, 3.8))
for r in range(6):
    s = [k0]
    while 0 < s[-1] < Ncap:
        s.append(s[-1] + rng.choice([-1, 1]))
    c = PALETTE[2] if s[-1] == Ncap else PALETTE[1]
    ax.step(np.arange(len(s)), s, where="post", color=c, lw=1.5, alpha=0.85)
    ax.plot(len(s) - 1, s[-1], "o", color=c, ms=6)
ax.axhline(0, color=PALETTE[1], lw=2.4)
ax.axhline(Ncap, color=PALETTE[2], lw=2.4)
ax.text(0.5, Ncap + 0.35, f"win: reach {Ncap}", color=PALETTE[2], fontsize=9.5)
ax.text(0.5, 0.35, "ruin: reach 0", color=PALETTE[1], fontsize=9.5)
ax.set_xlabel("bet number")
ax.set_ylabel("fortune")
ax.set_ylim(-0.6, Ncap + 1.2)
ax.set_title(f"Fair coin, start at {k0}: P(win) = {k0}/{Ncap}, mean duration {k0}·{Ncap-k0} = {k0*(Ncap-k0)}")
save(fig, "ch34-gamblers-ruin")
