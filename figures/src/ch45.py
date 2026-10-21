"""Chapter 45 figures: hawk-dove replicator dynamics and rock-paper-scissors cycles on the simplex."""
import sys, pathlib; sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
from style import plt, np, PALETTE, INK, MUTED, new_figure, diagram, node, box, arrow, save
from scipy.integrate import solve_ivp

# ---------------------------------------------------------------- hawk-dove
V, C = 4.0, 10.0
pstar = V / C
fig, (a1, a2) = new_figure(1, 2, figsize=(9.0, 3.8))
p = np.linspace(0, 1, 400)
pdot = p * (1 - p) * (V - C * p) / 2
a1.axhline(0, color=INK, lw=0.8)
a1.plot(p, pdot, color=PALETTE[0])
a1.fill_between(p, pdot, 0, where=pdot > 0, color=PALETTE[2], alpha=0.12)
a1.fill_between(p, pdot, 0, where=pdot < 0, color=PALETTE[1], alpha=0.12)
for q, fc in [(0, "white"), (1, "white"), (pstar, PALETTE[0])]:
    a1.plot(q, 0, "o", ms=8, mfc=fc, mec=PALETTE[0], mew=2, zorder=5)
for x0, d in [(0.12, 1), (0.25, 1), (0.6, -1), (0.8, -1)]:
    a1.annotate("", xy=(x0 + 0.06 * d, -0.045), xytext=(x0, -0.045),
                arrowprops=dict(arrowstyle="-|>", color=MUTED, lw=1.6))
a1.text(pstar + 0.03, 0.02, r"$p^*=V/C=0.4$", ha="left", color=INK)
a1.text(0.2, 0.06, "hawks spread", color=PALETTE[2], ha="center")
a1.text(0.72, -0.13, "hawks decline", color=PALETTE[1], ha="center")
a1.set_xlabel("share of hawks $p$")
a1.set_ylabel(r"$\dot p$")
a1.set_title("Phase line: one attractor inside")
a1.set_xlim(-0.03, 1.03)

t = np.linspace(0, 12, 400)
for k, p0 in enumerate([0.02, 0.15, 0.3, 0.55, 0.75, 0.97]):
    s = solve_ivp(lambda _, y: y * (1 - y) * (V - C * y) / 2, (0, 12), [p0], t_eval=t, rtol=1e-9)
    a2.plot(s.t, s.y[0], color=PALETTE[k % 6], lw=1.8)
a2.axhline(pstar, color=INK, ls="--", lw=1)
a2.text(12, pstar + 0.03, "ESS  $p^*=0.4$", ha="right", color=INK)
a2.set_xlabel("time $t$ (generations)")
a2.set_ylabel("share of hawks $p(t)$")
a2.set_title("Every mixed start converges")
a2.set_ylim(0, 1)
save(fig, "ch45-hawk-dove")

# ---------------------------------------------------------------- RPS on the simplex
def tern(x):
    x = np.atleast_2d(x)
    return x[:, 1] + 0.5 * x[:, 2], (np.sqrt(3) / 2) * x[:, 2]

def rps(a, b):
    A = np.array([[0, -b, a], [a, 0, -b], [-b, a, 0]], float)
    def f(_, x):
        fit = A @ x
        return x * (fit - x @ fit)
    return f

fig, axes = new_figure(1, 3, figsize=(10.5, 3.9))
cases = [(1, 1, "win = loss: closed cycles", [0.08, 0.16, 0.24]),
         (1.5, 1, "win > loss: spiral in", [0.28]),
         (1, 1.5, "win < loss: spiral out", [0.03])]
labels = ["Rock", "Paper", "Scissors"]
corners = np.eye(3)
for ax, (a, b, title, offs) in zip(axes, cases):
    cx, cy = tern(corners)
    ax.fill(np.r_[cx, cx[0]], np.r_[cy, cy[0]], color="#f8fafc", ec=INK, lw=1.2)
    for (u, v), lab, dx, dy in zip(zip(cx, cy), labels, [-0.02, 0.02, 0], [-0.06, -0.06, 0.04]):
        ax.text(u + dx, v + dy, lab, ha=["right", "left", "center"][labels.index(lab)], va="center", fontsize=9.5)
    T = 60 if a == b else (40 if a > b else 70)
    for k, d in enumerate(offs):
        x0 = np.array([1 / 3 + d, 1 / 3 - d / 2, 1 / 3 - d / 2])
        s = solve_ivp(rps(a, b), (0, T), x0, rtol=1e-10, atol=1e-12, max_step=0.05)
        u, v = tern(s.y.T)
        ax.plot(u, v, color=PALETTE[k if a == b else (2 if a > b else 1)], lw=1.4)
        i = len(u) // 7
        ax.annotate("", xy=(u[i + 3], v[i + 3]), xytext=(u[i], v[i]),
                    arrowprops=dict(arrowstyle="-|>", color=INK, lw=1.2))
    u, v = tern(np.array([1 / 3, 1 / 3, 1 / 3]))
    ax.plot(u, v, "o", ms=6, color=INK)
    ax.set_title(title, fontsize=10.5)
    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_xlim(-0.2, 1.2)
    ax.set_ylim(-0.12, 0.95)
save(fig, "ch45-rps-simplex")
