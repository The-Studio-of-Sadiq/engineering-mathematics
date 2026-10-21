"""Figures for Chapter 47 — Convexity: The Geometry of the Best Answer."""
import sys, pathlib; sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
from style import plt, np, PALETTE, INK, MUTED, new_figure, diagram, node, box, arrow, save

BLUE, RED, GREEN, PURPLE, ORANGE = PALETTE[:5]

# ---------------------------------------------------------------- figure 1: sets and separation
fig, axs = new_figure(1, 3, figsize=(10.5, 3.7))
for ax in axs:
    ax.set_aspect("equal"); ax.axis("off")

# (a) convex: an ellipse with a chord inside
ax = axs[0]
t = np.linspace(0, 2 * np.pi, 300)
ex, ey = 2.0 * np.cos(t), 1.25 * np.sin(t)
c, s = np.cos(0.4), np.sin(0.4)
X, Y = c * ex - s * ey, s * ex + c * ey
ax.fill(X, Y, color=BLUE, alpha=0.15); ax.plot(X, Y, color=BLUE)
p, q = np.array([-1.3, -0.9]), np.array([1.2, 1.1])
ax.plot(*zip(p, q), color=INK, lw=1.6); ax.plot(*zip(p, q), "o", color=INK, ms=5)
ax.text(p[0] - 0.15, p[1] - 0.3, r"$\mathbf{x}$", ha="center", color=INK)
ax.text(q[0] + 0.2, q[1] + 0.15, r"$\mathbf{y}$", ha="center", color=INK)
ax.set_title("convex: every chord stays inside", fontsize=10.5)
ax.set_xlim(-2.6, 2.6); ax.set_ylim(-2.2, 2.2)

# (b) non-convex: a crescent; a chord leaves the set
ax = axs[1]
th = np.linspace(0, 2 * np.pi, 400)
from matplotlib.path import Path as MPath
outer = np.c_[1.8 * np.cos(th), 1.8 * np.sin(th)]
inner = np.c_[0.9 + 1.4 * np.cos(th), 0.3 + 1.4 * np.sin(th)]
gx, gy = np.meshgrid(np.linspace(-2, 2, 400), np.linspace(-2, 2, 400))
inside = (gx**2 + gy**2 <= 1.8**2) & ((gx - 0.9)**2 + (gy - 0.3)**2 > 1.4**2)
ax.contourf(gx, gy, inside.astype(float), levels=[0.5, 1.5], colors=[RED], alpha=0.15)
ax.contour(gx, gy, inside.astype(float), levels=[0.5], colors=[RED], linewidths=2)
p, q = np.array([0.2, 1.6]), np.array([0.1, -1.65])
ax.plot(*zip(p, q), color=INK, lw=1.6, ls="--"); ax.plot(*zip(p, q), "o", color=INK, ms=5)
ax.text(0.45, 0.0, "chord leaves\nthe set", fontsize=9, color=MUTED, ha="left", va="center")
ax.set_title("non-convex: a chord escapes", fontsize=10.5)
ax.set_xlim(-2.6, 2.6); ax.set_ylim(-2.2, 2.2)

# (c) two disjoint convex sets, closest pair, separating hyperplane
ax = axs[2]
A = np.c_[-1.2 + 0.9 * np.cos(th), 0.6 + 0.9 * np.sin(th)]
B = np.array([[0.7, -0.6], [2.1, -0.3], [2.3, -1.8], [1.0, -1.9], [0.7, -0.6]])
ax.fill(A[:, 0], A[:, 1], color=BLUE, alpha=0.15); ax.plot(A[:, 0], A[:, 1], color=BLUE)
ax.fill(B[:, 0], B[:, 1], color=GREEN, alpha=0.15); ax.plot(B[:, 0], B[:, 1], color=GREEN)
b0 = B[0]; ca = np.array([-1.2, 0.6]); d = (b0 - ca) / np.linalg.norm(b0 - ca)
a0 = ca + 0.9 * d
ax.plot(*zip(a0, b0), "o-", color=INK, ms=4, lw=1.3)
mid = (a0 + b0) / 2; nrm = b0 - a0; tang = np.array([-nrm[1], nrm[0]]); tang /= np.linalg.norm(tang)
L = np.array([mid - 2.4 * tang, mid + 2.4 * tang])
ax.plot(L[:, 0], L[:, 1], color=PURPLE, lw=2)
base = mid + 1.3 * tang
ax.annotate("", xy=base + 0.8 * nrm / np.linalg.norm(nrm), xytext=base,
            arrowprops=dict(arrowstyle="-|>", color=PURPLE, lw=1.6))
ax.text(*(base + 0.95 * nrm / np.linalg.norm(nrm)), r"$\mathbf{a}$", color=PURPLE, fontsize=11)
ax.text(mid[0] - 0.4, mid[1] - 0.45, "closest pair", fontsize=8.5, color=MUTED, ha="right")
ax.text(-1.2, 0.6, r"$C$", ha="center", va="center", fontsize=12, color=BLUE)
ax.text(1.6, -1.2, r"$D$", ha="center", va="center", fontsize=12, color=GREEN)
ax.text(L[1, 0] + 0.05, L[1, 1] + 0.05, r"$\mathbf{a}^\mathsf{T}\mathbf{x}=b$", color=PURPLE, fontsize=10)
ax.set_title("disjoint convex sets: a hyperplane splits them", fontsize=10.5)
ax.set_xlim(-2.6, 3.0); ax.set_ylim(-2.4, 2.2)
save(fig, "ch47-sets-and-separation")

# ---------------------------------------------------------------- figure 2: tangents below a convex function
fig, axs = new_figure(1, 2, figsize=(10, 3.9))
x = np.linspace(-2.5, 3, 400)
f = lambda x: 0.3 * x**2 + np.exp(-0.8 * x)
df = lambda x: 0.6 * x - 0.8 * np.exp(-0.8 * x)
ax = axs[0]
ax.fill_between(x, f(x), 9, color=BLUE, alpha=0.10)
ax.text(-1.9, 7.6, "epigraph\n(convex set)", color=BLUE, fontsize=9.5)
for i, x0 in enumerate([-1.8, -0.6, 0.6, 1.8, 2.7]):
    ax.plot(x, f(x0) + df(x0) * (x - x0), color=MUTED, lw=1.1, ls="--")
    ax.plot(x0, f(x0), "o", color=ORANGE, ms=5, zorder=4)
ax.plot(x, f(x), color=BLUE, lw=2.4, label=r"$f(x)=0.3x^2+e^{-0.8x}$")
ax.set_ylim(-1, 9); ax.set_xlim(-2.5, 3)
ax.set_xlabel("$x$"); ax.set_title("convex: every tangent is a global under-estimator")
ax.legend(loc="upper center", fontsize=9)

g = lambda x: 0.5 * x**4 - 1.5 * x**2 + 0.3 * x + 1
dg = lambda x: 2 * x**3 - 3 * x + 0.3
ax = axs[1]
x = np.linspace(-1.9, 1.9, 400)
x0 = 0.25
ax.plot(x, g(x0) + dg(x0) * (x - x0), color=RED, lw=1.3, ls="--", label="tangent at $x_0$ cuts the graph")
ax.plot(x, g(x), color=INK, lw=2.4, label=r"$g(x)=0.5x^4-1.5x^2+0.3x+1$")
ax.plot(x0, g(x0), "o", color=RED, ms=6, zorder=4)
xs = np.linspace(-1.9, 1.9, 400); below = g(xs) < g(x0) + dg(x0) * (xs - x0)
ax.fill_between(xs, g(xs), g(x0) + dg(x0) * (xs - x0), where=below, color=RED, alpha=0.15)
ax.set_ylim(-1, 4); ax.set_xlabel("$x$")
ax.set_title("non-convex: the tangent lies above the graph")
ax.legend(loc="upper center", fontsize=9)
save(fig, "ch47-tangents")
