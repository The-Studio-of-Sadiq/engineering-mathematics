"""Figures for Chapter 51 — Optimal Control: Choosing the Best Future."""
import sys, pathlib; sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
from style import plt, np, PALETTE, INK, MUTED, new_figure, diagram, node, box, arrow, save

BLUE, RED, GREEN, PURPLE, ORANGE, TEAL, GREY = PALETTE

# ---------------------------------------------------------------- figure 1: bang-bang vs smooth
d, amax = 100.0, 1.0                        # 100 m, 1 m/s^2
Tb = 2 * np.sqrt(d / amax)                  # 20 s
Ts = np.sqrt(6 * d / amax)                  # 24.49 s, slowest-but-feasible smooth plan


def bang(t):
    ts = Tb / 2
    u = np.where(t < ts, amax, -amax)
    v = np.where(t < ts, amax * t, amax * (Tb - t))
    x = np.where(t < ts, 0.5 * amax * t**2, d - 0.5 * amax * (Tb - t)**2)
    return u, v, x


def smooth(t, T):
    c = 6 * d / T**2
    u = c * (1 - 2 * t / T)
    v = c * (t - t**2 / T)
    x = c * (t**2 / 2 - t**3 / (3 * T))
    return u, v, x


fig, axs = new_figure(2, 2, figsize=(10, 6.2))
tb = np.linspace(0, Tb, 801)
t2 = np.linspace(0, Tb, 400)
ts = np.linspace(0, Ts, 400)
ub, vb, xb = bang(tb)
u2, v2, x2 = smooth(t2, Tb)
us, vs, xs = smooth(ts, Ts)
specs = [(axs[0, 0], 0, "control u(t)  (m/s²)"), (axs[0, 1], 1, "speed v(t)  (m/s)"),
         (axs[1, 0], 2, "position x(t)  (m)")]
for ax, k, title in specs:
    ax.plot(tb, [ub, vb, xb][k], color=RED, lw=2.4, label="bang-bang, T = 20 s")
    ax.plot(t2, [u2, v2, x2][k], color=MUTED, lw=1.6, ls="--", label="min-energy, T = 20 s")
    ax.plot(ts, [us, vs, xs][k], color=BLUE, lw=2.0, label="min-energy, T = 24.5 s")
    ax.set_title(title, fontsize=10.5, loc="left")
    ax.set_xlim(0, 26)
axs[0, 0].axhspan(1, 1.7, color=RED, alpha=0.07); axs[0, 0].axhspan(-1.7, -1, color=RED, alpha=0.07)
axs[0, 0].text(13.5, 1.35, "forbidden: |u| > 1", fontsize=8.5, color=RED, ha="left", va="center")
axs[0, 0].set_ylim(-1.7, 1.7)
axs[0, 0].legend(loc="lower left", fontsize=8.2)
axs[1, 0].set_xlabel("time t (s)")

ax = axs[1, 1]
lam = -1 + tb / (Tb / 2)                    # costate lambda_2 (scaled), zero at the switch
ax.plot(tb, lam, color=PURPLE, lw=2.2, label=r"costate $\lambda_2(t)$")
ax.plot(tb, -np.sign(lam + 1e-12), color=RED, lw=1.4, ls=":", label=r"$u=-\mathrm{sign}\,\lambda_2$")
ax.axhline(0, color=INK, lw=0.8)
ax.plot(Tb / 2, 0, "o", color=INK, ms=6, zorder=5)
ax.annotate("switch at t = 10 s", xy=(Tb / 2, 0), xytext=(12.5, -0.65), fontsize=9,
            arrowprops=dict(arrowstyle="-|>", color=MUTED))
ax.set_title(r"why it switches: the costate crosses zero", fontsize=10.5, loc="left")
ax.set_xlim(0, 26); ax.set_ylim(-1.4, 1.4); ax.set_xlabel("time t (s)")
ax.legend(loc="lower right", fontsize=8.5)
save(fig, "ch51-bang-vs-smooth")

# ---------------------------------------------------------------- figure 2: phase plane and switching curve
fig, ax = new_figure(figsize=(7.4, 5.0))
v = np.linspace(-4, 4, 400)
ax.plot(-0.5 * v * np.abs(v), v, color=INK, lw=2.4, label=r"switching curve $x=-\frac{1}{2}v\,|v|$")
# background parabolas of u = +1 (x = v^2/2 + c) and u = -1 (x = -v^2/2 + c)
for c in np.arange(-10, 11, 2):
    ax.plot(0.5 * v**2 + c, v, color=BLUE, lw=0.6, alpha=0.25)
    ax.plot(-0.5 * v**2 + c, v, color=RED, lw=0.6, alpha=0.25)


def trajectory(x0, v0, dt=0.002):
    xs, vs = [x0], [v0]
    x, vv = x0, v0
    for _ in range(20000):
        s = x + 0.5 * vv * abs(vv)
        u = -np.sign(s) if abs(s) > 1e-9 else -np.sign(vv)
        if abs(s) < 5e-3:                   # on the curve: ride it home
            u = -np.sign(vv)
        vv += u * dt; x += vv * dt
        xs.append(x); vs.append(vv)
        if x * x + vv * vv < 1e-4:
            break
    return np.array(xs), np.array(vs)


starts = [(-6, 0), (-1.5, 3.4), (5, 1), (2, -3.2), (6, -1.5)]
for (x0, v0), col in zip(starts, [GREEN, ORANGE, PURPLE, TEAL, BLUE]):
    X, V = trajectory(x0, v0)
    ax.plot(X, V, color=col, lw=2)
    ax.plot(x0, v0, "o", color=col, ms=6, zorder=5)
    mid = len(X) // 4
    ax.annotate("", xy=(X[mid + 5], V[mid + 5]), xytext=(X[mid], V[mid]),
                arrowprops=dict(arrowstyle="-|>", color=col, lw=1.8))
ax.plot(0, 0, "*", color=INK, ms=14, zorder=6)
ax.text(0.25, -0.45, "target", fontsize=9, color=INK)
ax.text(-7.7, -3.7, "left of the curve:\nu = +1 (full throttle)", fontsize=9, color=BLUE)
ax.text(3.2, 2.2, "right of the curve:\nu = -1 (full brake)", fontsize=9, color=RED)
ax.set_xlim(-8, 8); ax.set_ylim(-4, 4)
ax.set_xlabel("position error x"); ax.set_ylabel("velocity v")
ax.legend(loc="upper right", fontsize=9)
save(fig, "ch51-phase-plane")
