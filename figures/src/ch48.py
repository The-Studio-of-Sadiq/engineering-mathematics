"""Figures for Chapter 48 — Duality: Solving Problems from the Other Side."""
import sys, pathlib; sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
from style import plt, np, PALETTE, INK, MUTED, new_figure, diagram, node, box, arrow, save

BLUE, RED, GREEN, PURPLE, ORANGE = PALETTE[:5]

# ------------------------------------------------- figure 1: primal and dual feasible regions
fig, axs = new_figure(1, 2, figsize=(10, 4.2))

ax = axs[0]
poly = np.array([[0, 0], [4, 0], [2, 3], [0, 4]])
ax.fill(poly[:, 0], poly[:, 1], color=BLUE, alpha=0.15)
ax.plot(*np.r_[poly, poly[:1]].T, color=BLUE, lw=2)
xs = np.linspace(-0.3, 5, 2)
for z, ls, col in [(12, ":", MUTED), (21, "--", RED)]:
    ax.plot(xs, (z - 3 * xs) / 5, ls=ls, color=col, lw=1.4)
ax.text(0.1, 2.45, r"$3x_1+5x_2=12$", color=MUTED, fontsize=9)
ax.text(2.95, 2.75, r"$3x_1+5x_2=21$", color=RED, fontsize=9)
ax.plot(2, 3, "o", color=RED, ms=7, zorder=5)
ax.annotate(r"$\mathbf{x}^*=(2,3)$, value 21", xy=(2, 3), xytext=(2.7, 4.2), fontsize=9.5,
            arrowprops=dict(arrowstyle="-", color=INK, lw=0.8))
ax.text(0.9, 1.1, "primal\nfeasible", color=BLUE, fontsize=10, ha="center")
ax.text(3.8, 0.4, r"$3x_1+2x_2\leq12$", fontsize=8.5, color=INK, rotation=-56)
ax.text(0.15, 3.2, r"$x_1+2x_2\leq8$", fontsize=8.5, color=INK, rotation=-27)
ax.set_xlim(-0.3, 5); ax.set_ylim(-0.3, 5); ax.set_aspect("equal")
ax.set_xlabel("$x_1$"); ax.set_ylabel("$x_2$")
ax.set_title("primal: maximise $3x_1+5x_2$", fontsize=11)

ax = axs[1]
Y1, Y2 = np.meshgrid(np.linspace(0, 4, 400), np.linspace(0, 2, 400))
feas = (Y1 + 3 * Y2 >= 3) & (2 * Y1 + 2 * Y2 >= 5)
ax.contourf(Y1, Y2, feas.astype(float), levels=[0.5, 1.5], colors=[GREEN], alpha=0.15)
y1 = np.linspace(0, 4, 200)
ax.plot(y1, (3 - y1) / 3, color=GREEN, lw=1.3)
ax.plot(y1, (5 - 2 * y1) / 2, color=GREEN, lw=1.3)
for z, ls, col in [(30, ":", MUTED), (21, "--", RED)]:
    ax.plot(y1, (z - 8 * y1) / 12, ls=ls, color=col, lw=1.4)
ax.text(3.05, 0.62, r"$8y_1+12y_2=30$", color=MUTED, fontsize=9)
ax.text(0.15, 1.45, r"$8y_1+12y_2=21$", color=RED, fontsize=9)
ax.plot(2.25, 0.25, "o", color=RED, ms=7, zorder=5)
ax.annotate(r"$\mathbf{y}^*=(2.25,\,0.25)$, value 21", xy=(2.25, 0.25), xytext=(2.3, 1.45), fontsize=9.5,
            arrowprops=dict(arrowstyle="-", color=INK, lw=0.8))
ax.text(3.3, 1.7, "dual\nfeasible", color=GREEN, fontsize=10, ha="center")
ax.set_xlim(0, 4); ax.set_ylim(0, 2)
ax.set_xlabel("$y_1$ (price of resource 1)"); ax.set_ylabel("$y_2$ (price of resource 2)")
ax.set_title("dual: minimise $8y_1+12y_2$", fontsize=11)
save(fig, "ch48-primal-dual")

# ------------------------------------------------- figure 2: the image set and the duality gap
fig, axs = new_figure(1, 2, figsize=(10, 4.2))

# convex: min x^2 s.t. 1 - x <= 0 ; u = 1 - x, t = x^2 = (1-u)^2
ax = axs[0]
u = np.linspace(-1.6, 1.0, 300)
ax.fill_between(u, (1 - u) ** 2, 8, color=BLUE, alpha=0.10)
ax.plot(u, (1 - u) ** 2, color=BLUE, lw=2.2, label=r"$\mathcal{G}=\{(g(x),f(x))\}$")
ax.plot(u, 1 - 2 * u, color=PURPLE, lw=1.6, ls="--", label=r"supporting line, slope $-\lambda^*=-2$")
ax.axvline(0, color=MUTED, lw=1)
ax.plot(0, 1, "o", color=RED, ms=7, zorder=5)
ax.text(0.08, 1.25, r"$p^*=d^*=1$", color=RED, fontsize=10)
ax.text(-1.5, 0.4, r"feasible side $u\leq0$", fontsize=9, color=MUTED)
ax.set_xlim(-1.6, 1.0); ax.set_ylim(-0.5, 6.5)
ax.set_xlabel(r"$u=g(x)$ (constraint value)"); ax.set_ylabel(r"$t=f(x)$ (cost)")
ax.set_title("convex problem: no gap", fontsize=11)
ax.legend(loc="upper right", fontsize=8.5)

# non-convex: min -x^2 s.t. x - 1 <= 0, x in [0,2] ; u = x - 1, t = -(u+1)^2
ax = axs[1]
u = np.linspace(-1, 1, 300)
ax.plot(u, -(u + 1) ** 2, color=INK, lw=2.2, label=r"$\mathcal{G}$ (a curved arc)")
ax.plot(np.linspace(-1.3, 1.3, 2), -2 - 2 * np.linspace(-1.3, 1.3, 2), color=PURPLE, lw=1.6, ls="--",
        label=r"best supporting line, slope $-2$")
ax.axvline(0, color=MUTED, lw=1)
ax.plot(0, -1, "o", color=RED, ms=7, zorder=5); ax.text(0.07, -0.85, r"$p^*=-1$", color=RED, fontsize=10)
ax.plot(0, -2, "s", color=PURPLE, ms=7, zorder=5); ax.text(0.1, -1.85, r"$d^*=-2$", color=PURPLE, fontsize=10)
ax.annotate("", xy=(-0.04, -1.05), xytext=(-0.04, -1.95),
            arrowprops=dict(arrowstyle="<->", color=ORANGE, lw=1.6))
ax.text(-0.55, -1.55, "duality\ngap = 1", color=ORANGE, fontsize=9.5, ha="center")
ax.plot([-1, 1], [0, -4], "o", color=INK, ms=5)
ax.set_xlim(-1.3, 1.3); ax.set_ylim(-4.6, 0.8)
ax.set_xlabel(r"$u=g(x)$"); ax.set_ylabel(r"$t=f(x)$")
ax.set_title("non-convex problem: a gap opens", fontsize=11)
ax.legend(loc="lower left", fontsize=8.5)
save(fig, "ch48-duality-gap")
