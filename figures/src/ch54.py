"""Figures for Chapter 54 — Robust and Stochastic Optimization."""
import sys, pathlib; sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
from style import plt, np, PALETTE, INK, MUTED, new_figure, diagram, node, box, arrow, save
from scipy.optimize import minimize

BLUE, RED, GREEN, PURPLE, ORANGE, TEAL, GREY = PALETTE

# Running example: maximise x1 + x2  s.t.  2x1 + x2 <= 8,  x1 + 3x2 <= 9,  x >= 0,
# with every coefficient of A perturbed by independent N(0, sigma^2) noise.
A = np.array([[2.0, 1.0], [1.0, 3.0]])
b = np.array([8.0, 9.0])
c = np.array([1.0, 1.0])
SIG = 0.1


def robust(rho):
    cons = [{"type": "ineq", "fun": lambda x, i=i: b[i] - A[i] @ x - rho * np.linalg.norm(x)} for i in range(2)]
    r = minimize(lambda x: -c @ x, [1, 1], constraints=cons, bounds=[(0, None)] * 2,
                 method="SLSQP", options={"ftol": 1e-12})
    return r.x


def region_mask(X, Y, rho):
    R = np.hypot(X, Y)
    return (2 * X + Y + rho * R <= 8) & (X + 3 * Y + rho * R <= 9) & (X >= 0) & (Y >= 0)


# ------------------------------------------------- figure 1: scenarios and shrinking feasible regions
fig, axs = new_figure(1, 2, figsize=(10.5, 4.6))
xs = np.linspace(0, 4.6, 400)

ax = axs[0]
rng = np.random.default_rng(54)
for k in range(40):
    Ap = A + SIG * rng.standard_normal((2, 2))
    ax.plot(xs, (b[0] - Ap[0, 0] * xs) / Ap[0, 1], color=BLUE, lw=0.6, alpha=0.35)
    ax.plot(xs, (b[1] - Ap[1, 0] * xs) / Ap[1, 1], color=ORANGE, lw=0.6, alpha=0.35)
poly = np.array([[0, 0], [4, 0], [3, 2], [0, 3]])
ax.fill(poly[:, 0], poly[:, 1], color=GREY, alpha=0.10)
ax.plot(*np.r_[poly, poly[:1]].T, color=INK, lw=1.8)
ax.plot(3, 2, "o", color=RED, ms=8, zorder=5)
ax.annotate("nominal optimum $(3,2)$\nviolates a sampled\nconstraint 75% of the time", xy=(3, 2),
            xytext=(3.25, 3.25), fontsize=9, arrowprops=dict(arrowstyle="-", color=INK, lw=0.8))
ax.text(0.9, 1.0, "nominal\nfeasible set", ha="center", color=INK, fontsize=9.5)
ax.text(0.12, 4.0, "blue: samples of $2x_1+x_2\\leq8$", color=BLUE, fontsize=9)
ax.text(0.12, 3.75, "orange: samples of $x_1+3x_2\\leq9$", color=ORANGE, fontsize=9)
ax.set_xlim(0, 4.6); ax.set_ylim(0, 4.2); ax.set_aspect("equal")
ax.set_xlabel("$x_1$"); ax.set_ylabel("$x_2$")
ax.set_title("40 sampled versions of the constraints ($\\sigma=0.1$)", fontsize=11)

ax = axs[1]
X, Y = np.meshgrid(np.linspace(0, 4.6, 600), np.linspace(0, 4.2, 600))
rhos = [0.0, 0.15, 0.3, 0.45]
cols = [GREY, BLUE, PURPLE, GREEN]
for rho, col in zip(rhos, cols):
    m = region_mask(X, Y, rho).astype(float)
    ax.contourf(X, Y, m, levels=[0.5, 1.5], colors=[col], alpha=0.10)
    ax.contour(X, Y, m, levels=[0.5], colors=[col], linewidths=1.6)
path = np.array([robust(r) for r in np.linspace(0, 0.45, 30)])
ax.plot(path[:, 0], path[:, 1], color=RED, lw=1.4, ls="--")
for rho, col in zip(rhos, cols):
    x = robust(rho)
    ax.plot(*x, "o", color=col if rho else RED, ms=7, zorder=5, mec="white",
            label=f"$\\rho={rho:g}$: optimum $({x[0]:.2f},\\,{x[1]:.2f})$, value {c @ x:.2f}")
ax.plot(xs, 5 - xs, color=MUTED, lw=1, ls=":")
ax.text(3.6, 1.6, "$x_1+x_2=5$", color=MUTED, fontsize=8.8, rotation=-45)
ax.legend(loc="upper right", fontsize=8.3, frameon=True, framealpha=0.95)
ax.set_xlim(0, 4.6); ax.set_ylim(0, 4.2); ax.set_aspect("equal")
ax.set_xlabel("$x_1$"); ax.set_ylabel("$x_2$")
ax.set_title(r"robust sets $\bar{\mathbf{a}}_i^\mathsf{T}\mathbf{x}+\rho\Vert\mathbf{x}\Vert\leq b_i$", fontsize=11)
save(fig, "ch54-region")

# ------------------------------------------------- figure 2: the price of robustness
M = 400_000
Ap = A[None] + SIG * rng.standard_normal((M, 2, 2))
rho_grid = np.linspace(0, 0.5, 26)
val, pv = [], []
for rho in rho_grid:
    x = robust(rho)
    val.append(c @ x)
    pv.append(max(((Ap @ x) > b).any(1).mean(), 1.0 / M))
val, pv = np.array(val), np.array(pv)
for r in [0, 0.1645, 0.3]:
    i = np.argmin(abs(rho_grid - r)) if r != 0.1645 else None
    x = robust(r)
    print(f"rho={r}: x={x.round(4)}, value={c @ x:.4f}, P(violate)={((Ap @ x) > b).any(1).mean():.4f}")

fig, ax = new_figure(figsize=(7.2, 4.2))
ax.semilogx(pv, val, color=BLUE, lw=2)
ax.plot(pv, val, "o", color=BLUE, ms=3.5)
marks = [(0.0, "nominal, $\\rho=0$"), (0.1645, "$\\rho=\\sigma z_{0.95}=0.164$"), (0.3, "$\\rho=0.3$")]
for r, lab in marks:
    x = robust(r)
    p = max(((Ap @ x) > b).any(1).mean(), 1.0 / M)
    ax.plot(p, c @ x, "o", color=RED, ms=8, zorder=5)
    ax.annotate(lab + f"\nvalue {c @ x:.2f}, P(fail) {p:.3f}", xy=(p, c @ x), xytext=(p * 1.6, c @ x - 0.25),
                fontsize=8.8, arrowprops=dict(arrowstyle="-", color=INK, lw=0.7))
ax.set_xlabel("probability that the plan violates a constraint (Monte Carlo)")
ax.set_ylabel("objective value $x_1+x_2$")
ax.set_title("the price of robustness", fontsize=11)
ax.set_xlim(1e-5, 1.5); ax.set_ylim(4.0, 5.1)
save(fig, "ch54-tradeoff")
