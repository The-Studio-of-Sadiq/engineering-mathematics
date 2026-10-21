"""Figures for Chapter 92 — Mathematical Physics: The Structure Behind Nature."""
import sys, pathlib; sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
from style import plt, np, PALETTE, INK, MUTED, new_figure, diagram, node, box, arrow, save

BLUE, RED, GREEN, PURPLE, ORANGE, TEAL, GREY = PALETTE
g = 9.81


# ---------------------------------------------------------------- figure 1: pendulum phase portrait
def rhs(s):
    q, p = s
    return np.array([p, -np.sin(q)])


def rk4_flow(s, t, n=3000):
    h = t / n
    for _ in range(n):
        k1 = rhs(s); k2 = rhs(s + h / 2 * k1); k3 = rhs(s + h / 2 * k2); k4 = rhs(s + h * k3)
        s = s + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
    return s


def shoelace(q, p):
    return 0.5 * abs(np.dot(q, np.roll(p, -1)) - np.dot(p, np.roll(q, -1)))


fig, ax = new_figure(figsize=(7.6, 4.6))
Q, P = np.meshgrid(np.linspace(-2 * np.pi, 2 * np.pi, 600), np.linspace(-3.2, 3.2, 400))
Hm = 0.5 * P**2 - np.cos(Q)
ax.contour(Q, P, Hm, levels=[-0.9, -0.6, -0.2, 0.3, 0.7], colors=BLUE, linewidths=1.2, linestyles="solid")
ax.contour(Q, P, Hm, levels=[1.5, 2.3, 3.5], colors=GREEN, linewidths=1.2)
ax.contour(Q, P, Hm, levels=[1.0], colors=RED, linewidths=2.2)
for qe in (-2 * np.pi, 0, 2 * np.pi):
    ax.plot(qe, 0, "o", color=INK, ms=5, zorder=5)
for qe in (-np.pi, np.pi):
    ax.plot(qe, 0, "o", mfc="white", mec=RED, mew=2, ms=7, zorder=5)

# Liouville: a disc of initial conditions carried by the exact flow
th = np.linspace(0, 2 * np.pi, 400, endpoint=False)
c0, r0 = np.array([-2.6, 0.9]), 0.32
blob0 = np.c_[c0[0] + r0 * np.cos(th), c0[1] + r0 * np.sin(th)].T
areas = []
for t, col in [(0.0, ORANGE), (2.0, PURPLE), (4.0, TEAL)]:
    b = blob0 if t == 0 else rk4_flow(blob0, t)
    areas.append(shoelace(*b))
    ax.fill(b[0], b[1], color=col, alpha=0.45, zorder=4)
    ax.plot(np.r_[b[0], b[0][0]], np.r_[b[1], b[1][0]], color=col, lw=1.2, zorder=4)
ax.text(-2.6, 1.42, "t = 0", color=ORANGE, ha="center", fontsize=9, fontweight="bold")
ax.annotate(f"same area at t = 0, 2, 4:  {areas[0]:.3f}, {areas[1]:.3f}, {areas[2]:.3f}",
            xy=(0.35, 2.05), xytext=(-1.0, 2.75), fontsize=9, color=INK,
            arrowprops=dict(arrowstyle="-|>", color=MUTED, lw=1.2))
ax.text(4.6, 2.55, "rotation (H > 1)", color=GREEN, fontsize=9, ha="center",
        bbox=dict(fc="white", ec="none", pad=1))
ax.text(0, -0.25, "libration\n(H < 1)", color=BLUE, fontsize=9, ha="center", va="top",
        bbox=dict(fc="white", ec="none", pad=1))
ax.text(np.pi, -0.42, "unstable\nequilibrium", color=RED, fontsize=8.5, ha="center", va="top",
        bbox=dict(fc="white", ec="none", pad=1))
ax.text(-4.75, -2.05, "separatrix H = 1", color=RED, fontsize=9, ha="center",
        bbox=dict(fc="white", ec="none", pad=1))
ax.set_xlim(-2 * np.pi, 2 * np.pi)
ax.set_ylim(-3.2, 3.2)
ax.set_xticks(np.arange(-2, 3) * np.pi)
ax.set_xticklabels([r"$-2\pi$", r"$-\pi$", "0", r"$\pi$", r"$2\pi$"])
ax.set_xlabel(r"angle $q$ (rad)")
ax.set_ylabel(r"momentum $p$")
ax.set_title(r"Pendulum phase space: level sets of $H=\frac{1}{2}p^2-\cos q$")
save(fig, "ch92-pendulum-phase-space")
print("blob areas", areas)


# ---------------------------------------------------------------- figure 2: brachistochrone
def travel(x, y):
    v = np.sqrt(2 * g * np.maximum(-y, 0))
    L = np.hypot(np.diff(x), np.diff(y))
    return np.sum(2 * L / (v[:-1] + v[1:]))


fig, ax = new_figure(figsize=(7.4, 4.0))
B = (np.pi, -2.0)
xs = np.linspace(0, np.pi, 4001); ys = -2 * xs / np.pi
R = (np.pi**2 + 4) / (2 * np.pi)
ph = np.linspace(np.pi, 2 * np.pi + np.arctan2(-2, np.pi - R), 4001)
xc, yc = R + R * np.cos(ph), R * np.sin(ph)
tt = np.linspace(0, np.pi, 4001)
xk, yk = tt - np.sin(tt), -(1 - np.cos(tt))
curves = [(xs, ys, GREY, "straight line"), (xc, yc, ORANGE, "circular arc"),
          (xk, yk, BLUE, "cycloid")]
for x, y, col, lab in curves:
    ax.plot(x, y, color=col, lw=2.4, label=f"{lab}: {travel(x, y):.3f} s")
    print(lab, travel(x, y))
# a few beads at equal times along the cycloid and the line
for tl in (0.25, 0.5, 0.75):
    t_c = tl * np.sqrt(g)                   # cycloid (radius 1): parameter = t * sqrt(g / r)
    ax.plot(t_c - np.sin(t_c), -(1 - np.cos(t_c)), "o", color=BLUE, ms=6, zorder=5)
    s = 0.5 * g * (2 / np.hypot(np.pi, 2)) * tl**2   # distance slid along the line
    u = s / np.hypot(np.pi, 2)
    ax.plot(u * np.pi, -2 * u, "o", color=GREY, ms=6, zorder=5)
ax.plot(0, 0, "o", color=INK, ms=7); ax.text(0.08, 0.08, "A (release)", fontsize=9.5)
ax.plot(*B, "o", color=INK, ms=7); ax.text(B[0] - 0.05, B[1] - 0.22, "B", fontsize=9.5, ha="right")
ax.text(1.95, -0.12, "dots: bead positions at\nt = 0.25, 0.50, 0.75 s", fontsize=8.5, color=MUTED)
ax.set_aspect("equal")
ax.set_xlim(-0.15, 3.6); ax.set_ylim(-2.4, 0.4)
ax.set_xlabel("x (m)"); ax.set_ylabel("y (m)")
ax.legend(loc="upper right", fontsize=9, bbox_to_anchor=(1.0, 0.82))
ax.set_title("Brachistochrone: fastest frictionless slide from A to B (g = 9.81 m/s²)")
save(fig, "ch92-brachistochrone")
