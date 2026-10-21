"""Figures for Chapter 89 — Partial Differential Equations Beyond the Basics."""
import sys, pathlib; sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
from style import plt, np, PALETTE, INK, MUTED, new_figure, diagram, node, box, arrow, save

BLUE, RED, GREEN, PURPLE, ORANGE, TEAL, GREY = PALETTE

# ---------------------------------------------------------------- figure 1: a kinked weak solution, judged by test functions
a = 0.4
x = np.linspace(0, 1, 1001)
u = np.where(x < a, x * (1 - a), a * (1 - x))
du = np.where(x < a, 1 - a, -a)


def bump(x, c, w):
    s = (x - c) / w
    out = np.zeros_like(x)
    m = np.abs(s) < 1
    out[m] = np.exp(-1 / (1 - s[m] ** 2)) * np.e
    return out


fig, axs = new_figure(1, 2, figsize=(10.4, 3.9))
ax = axs[0]
ax.plot(x, u, color=BLUE, lw=2.6, label=r"$u(x)$: string under a point load at $x=0.4$")
for c, w, col in [(0.18, 0.12, GREEN), (0.45, 0.15, ORANGE), (0.78, 0.13, PURPLE)]:
    v = 0.12 * bump(x, c, w)
    ax.fill_between(x, 0, v, color=col, alpha=0.18)
    ax.plot(x, v, color=col, lw=1.5)
ax.annotate("kink: no $u''$ here", xy=(a, a * (1 - a)), xytext=(0.55, 0.27),
            arrowprops=dict(arrowstyle="-|>", color=RED, lw=1.4), color=RED, fontsize=10)
ax.text(0.78, 0.135, "smooth test\nfunctions $v$", ha="center", va="bottom", fontsize=9.5, color=MUTED)
ax.axhline(0, color=INK, lw=0.8)
ax.set_xlim(0, 1); ax.set_ylim(-0.02, 0.34)
ax.set_xlabel("$x$"); ax.set_ylabel("displacement")
ax.legend(loc="upper left", fontsize=9)
ax.set_title("The solution and the probes that judge it", fontsize=11)

ax = axs[1]
ax.plot(x[x < a], du[x < a], color=BLUE, lw=2.6)
ax.plot(x[x > a], du[x > a], color=BLUE, lw=2.6, label=r"$u'(x)$: a step of height $-1$")
ax.plot([a, a], [-a, 1 - a], color=BLUE, lw=1.2, ls=":")
ax.axhline(0, color=INK, lw=0.8)
ax.annotate("", xy=(a + 0.02, -a), xytext=(a + 0.02, 1 - a),
            arrowprops=dict(arrowstyle="<->", color=RED, lw=1.4))
ax.text(a + 0.045, 0.1, r"jump $=-1$" "\n" r"so $-u''=\delta_{0.4}$", color=RED, fontsize=10, va="center")
ax.text(0.5, -0.75, r"weak form: $\int_0^1 u'v'\,dx = v(0.4)$ for every $v$",
        ha="center", fontsize=10.5, color=INK,
        bbox=dict(boxstyle="round,pad=0.35", fc="#f3f4f6", ec=MUTED))
ax.set_xlim(0, 1); ax.set_ylim(-0.95, 0.85)
ax.set_xlabel("$x$"); ax.set_ylabel("slope")
ax.legend(loc="upper right", fontsize=9)
ax.set_title("One derivative exists, the second does not", fontsize=11)
save(fig, "ch89-weak-solution")

# ---------------------------------------------------------------- figure 2: hat-function basis and an FEM approximation
n = 6
nodes = np.linspace(0, 1, n + 1)
fig, axs = new_figure(2, 1, figsize=(8.4, 5.6), sharex=True)
ax = axs[0]
cols = [BLUE, RED, GREEN, PURPLE, ORANGE]
for i in range(1, n):
    phi = np.interp(x, nodes, np.eye(n + 1)[i])
    ax.plot(x, phi, color=cols[(i - 1) % 5], lw=2)
    ax.text(nodes[i], 1.06, rf"$\varphi_{i}$", ha="center", fontsize=10, color=cols[(i - 1) % 5])
ax.plot(nodes, 0 * nodes, "o", color=INK, ms=5, zorder=5)
ax.set_ylim(-0.08, 1.25); ax.set_ylabel("basis function")
ax.set_title(r"Hat functions: $\varphi_i=1$ at node $i$, $0$ at every other node", fontsize=11)

ax = axs[1]
# two-material bar: -(k u')' = 1 on (0,1), k=1 on (0,0.5), k=5 on (0.5,1), u(0)=u(1)=0
k1, k2, m = 1.0, 5.0, 0.5
# exact: k u' = C - x ; u continuous at m, u(1)=0
C = (m**2 / (2 * k1) + (1 - m**2) / (2 * k2)) / (m / k1 + (1 - m) / k2)
uex = np.where(x < m, (C * x - x**2 / 2) / k1,
               (C * m - m**2 / 2) / k1 + (C * (x - m) - (x**2 - m**2) / 2) / k2)
un = np.interp(nodes, x, uex)
ax.plot(x, uex, color=GREY, lw=4, alpha=0.35, label="exact solution (kink at the material interface)")
ax.plot(nodes, un, "-o", color=BLUE, lw=2, ms=5, label=r"FEM: $u_h=\sum_i c_i\varphi_i$")
for i in range(1, n):
    ax.plot(x, un[i] * np.interp(x, nodes, np.eye(n + 1)[i]), color=cols[(i - 1) % 5], lw=1, ls="--", alpha=0.8)
ax.axvline(m, color=RED, lw=1, ls=":")
ax.text(m + 0.01, 0.005, "  $k=1$ | $k=5$", color=RED, fontsize=9.5, ha="center", va="bottom")
ax.set_xlim(0, 1); ax.set_ylim(-0.003, 0.105)
ax.set_xlabel("$x$"); ax.set_ylabel("$u$")
ax.legend(loc="upper right", fontsize=9)
ax.set_title(r"Two-material bar $-(k u')'=1$: dashed lines are the scaled hats $c_i\varphi_i$", fontsize=11)
save(fig, "ch89-hat-basis")

# ---------------------------------------------------------------- figure 3: corner singularity on an L-shaped domain
N = 401
g = np.linspace(-1, 1, N)
X, Y = np.meshgrid(g, g)
R = np.hypot(X, Y)
TH = np.mod(np.arctan2(Y, X), 2 * np.pi)
inside = ~((X > 0) & (Y < 0))
U = np.where(inside, R ** (2 / 3) * np.sin(2 * TH / 3), np.nan)
G = np.where(inside, (2 / 3) * np.maximum(R, 1e-3) ** (-1 / 3), np.nan)

fig, axs = new_figure(1, 2, figsize=(10.0, 4.5))
ax = axs[0]
im = ax.pcolormesh(X, Y, U, cmap="viridis", shading="auto")
ax.contour(X, Y, U, levels=np.linspace(0.1, 1.4, 10), colors="white", linewidths=0.8)
ax.plot([0, 1, 1, -1, -1, 0, 0], [0, 0, 1, 1, -1, -1, 0], color=INK, lw=1.6)
ax.set_aspect("equal"); ax.grid(False)
ax.set_title(r"$u=r^{2/3}\sin(2\theta/3)$ solves $\nabla^2u=0$", fontsize=11)
fig.colorbar(im, ax=ax, shrink=0.85, label="$u$")
ax.annotate("re-entrant corner\n(interior angle $270^\\circ$)", xy=(0, 0), xytext=(0.15, -0.6),
            color=RED, fontsize=9.5, arrowprops=dict(arrowstyle="-|>", color=RED, lw=1.3))

ax = axs[1]
from matplotlib.colors import PowerNorm
im = ax.pcolormesh(X, Y, G, cmap="magma", shading="auto", norm=PowerNorm(gamma=0.5, vmin=0.6, vmax=4))
ax.contour(X, Y, G, levels=[0.8, 1.0, 1.5, 2.0, 3.0], colors="white", linewidths=0.8)
ax.plot([0, 1, 1, -1, -1, 0, 0], [0, 0, 1, 1, -1, -1, 0], color=INK, lw=1.6)
ax.set_aspect("equal"); ax.grid(False)
ax.set_title(r"$|\nabla u|=\frac{2}{3}r^{-1/3}$: infinite at the corner", fontsize=11)
fig.colorbar(im, ax=ax, shrink=0.85, label=r"$|\nabla u|$ (colour clipped at 4)")
save(fig, "ch89-corner-singularity")
