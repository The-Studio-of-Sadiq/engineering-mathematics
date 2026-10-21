"""Chapter 84 figures: Lorenz attractor and divergence of nearby orbits (own RK4), logistic-map cobwebs."""
import sys, pathlib; sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
from style import plt, np, PALETTE, INK, MUTED, new_figure, diagram, node, box, arrow, save

SIGMA, RHO, BETA = 10.0, 28.0, 8.0 / 3.0


def lorenz(u):
    x, y, z = u
    return np.array([SIGMA * (y - x), x * (RHO - z) - y, x * y - BETA * z])


def rk4(u, h, n):
    out = np.empty((n + 1, 3)); out[0] = u
    for i in range(n):
        k1 = lorenz(u); k2 = lorenz(u + h / 2 * k1); k3 = lorenz(u + h / 2 * k2); k4 = lorenz(u + h * k3)
        u = u + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
        out[i + 1] = u
    return out


def lorenz_figures():
    h = 0.005
    u0 = rk4(np.array([1.0, 1.0, 1.0]), h, 4000)[-1]          # settle onto the attractor first
    n = int(30 / h)
    a = rk4(u0, h, n)
    b = rk4(u0 + np.array([1e-9, 0, 0]), h, n)
    t = np.arange(n + 1) * h

    # projection with the two trajectories
    fig, ax = new_figure(figsize=(7.2, 5.2))
    bg = rk4(a[-1], h, 20000)
    ax.plot(bg[:, 0], bg[:, 2], color="#d1d5db", lw=0.5)
    i0, i1 = int(19 / h), int(28 / h)
    ax.plot(a[i0:i1, 0], a[i0:i1, 2], color=PALETTE[0], lw=1.1, label="orbit A, t = 19 to 28")
    ax.plot(b[i0:i1, 0], b[i0:i1, 2], color=PALETTE[1], lw=1.1, label="orbit B (x(0) shifted by 10⁻⁹), same times")
    c = np.sqrt(BETA * (RHO - 1))
    ax.plot([c, -c], [RHO - 1, RHO - 1], "X", color=INK, ms=9, mec="white", zorder=5)
    ax.text(c - 1.5, RHO - 3.5, "C₊", fontsize=10, color=INK)
    ax.text(-c - 1.5, RHO - 3.5, "C₋", fontsize=10, color=INK)
    ax.plot(a[i0, 0], a[i0, 2], "o", color=INK, ms=7, zorder=6)
    ax.annotate("t = 19: still together", (a[i0, 0], a[i0, 2]), (a[i0, 0] - 13, a[i0, 2] + 10), fontsize=9,
                arrowprops=dict(arrowstyle="-", color=MUTED), bbox=dict(fc="white", ec="none", alpha=0.9))
    for o, col in ((a, PALETTE[0]), (b, PALETTE[1])):
        ax.plot(o[i1, 0], o[i1, 2], "o", color=col, ms=8, mec="white", mew=1.5, zorder=6)
    ax.text(-19, 3, "dots: where A and B are at t = 28", fontsize=9, color=MUTED)
    ax.set_xlabel("x"); ax.set_ylabel("z")
    ax.set_title("Lorenz attractor (σ = 10, ρ = 28, β = 8/3), x–z projection")
    ax.legend(loc="upper left", fontsize=8.5, frameon=True, framealpha=0.9)
    ax.set_ylim(0, 56)
    save(fig, "ch84-lorenz-attractor")

    # log separation vs time
    d = np.linalg.norm(a - b, axis=1)
    fig, ax = new_figure(figsize=(7.2, 4.0))
    ax.semilogy(t, d, color=PALETTE[0], lw=1.3, label="separation ‖A(t) − B(t)‖")
    tt = np.linspace(0, 23, 50)
    ax.semilogy(tt, 1e-9 * np.exp(0.906 * tt), color=PALETTE[1], ls="--", lw=1.8,
                label=r"$10^{-9}\,e^{0.906\,t}$  (slope $\lambda_1 = 0.906$)")
    ax.axhline(40, color=MUTED, ls=":", lw=1.4)
    ax.text(0.5, 60, "size of the attractor: separation saturates", fontsize=9, color=MUTED)
    ax.set_xlabel("time t"); ax.set_ylabel("separation (log scale)")
    ax.set_ylim(1e-11, 1e3); ax.set_xlim(0, 30)
    ax.set_title("A 10⁻⁹ difference grows exponentially, then saturates")
    ax.legend(loc="lower right", fontsize=9)
    save(fig, "ch84-lorenz-divergence")


def cobweb(ax, r, x0, n, title):
    xs = np.linspace(0, 1, 400)
    ax.plot(xs, r * xs * (1 - xs), color=PALETTE[0], lw=2, label="f(x) = r x (1 − x)")
    ax.plot(xs, xs, color=MUTED, lw=1.2, ls="--", label="y = x")
    x = x0
    px, py = [x], [0]
    for _ in range(n):
        y = r * x * (1 - x)
        px += [x, y]; py += [y, y]
        x = y
    ax.plot(px, py, color=PALETTE[1], lw=0.7)
    ax.plot(x0, 0, "o", color=PALETTE[1], ms=6, clip_on=False, zorder=5)
    ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.set_aspect("equal")
    ax.set_xlabel("$x_n$"); ax.set_ylabel("$x_{n+1}$")
    ax.set_title(title)


def cobwebs():
    fig, axes = new_figure(1, 3, figsize=(12, 4.3))
    cobweb(axes[0], 2.8, 0.08, 40, "r = 2.8: spiral into x* = 0.643")
    cobweb(axes[1], 3.3, 0.08, 60, "r = 3.3: settles on a 2-cycle")
    cobweb(axes[2], 3.9, 0.08, 50, "r = 3.9: never settles")
    axes[0].legend(loc="upper left", fontsize=8.5)
    save(fig, "ch84-logistic-cobweb")


if __name__ == "__main__":
    lorenz_figures()
    cobwebs()
