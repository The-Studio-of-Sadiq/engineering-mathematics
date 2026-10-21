"""Chapter 83 figures: nonlinear phase portraits (pendulum, predator-prey, Van der Pol limit cycle)."""
import sys, pathlib; sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
from style import plt, np, PALETTE, INK, MUTED, new_figure, diagram, node, box, arrow, save
from scipy.integrate import solve_ivp


def pendulum_portrait():
    fig, axes = new_figure(1, 2, figsize=(11, 4.4))
    th = np.linspace(-2 * np.pi - 0.6, 2 * np.pi + 0.6, 300)
    w = np.linspace(-3.2, 3.2, 240)
    T, W = np.meshgrid(th, w)
    for ax, b, title in [(axes[0], 0.0, "undamped:  θ'' + sin θ = 0"),
                         (axes[1], 0.25, "damped:  θ'' + 0.25 θ' + sin θ = 0")]:
        U, V = W, -np.sin(T) - b * W
        speed = np.hypot(U, V)
        ax.streamplot(T, W, U, V, color=speed, cmap="Blues", density=1.4, linewidth=0.9, arrowsize=0.9)
        if b == 0:
            # separatrix: E = w^2/2 + 1 - cos(theta) = 2
            ws = 2 * np.cos(th / 2)
            ax.plot(th, ws, color=PALETTE[1], lw=2.2, label="separatrix, E = 2")
            ax.plot(th, -ws, color=PALETTE[1], lw=2.2)
            ax.legend(loc="upper right", fontsize=9, frameon=True, framealpha=0.9)
        for k in (-2, 0, 2):
            ax.plot(k * np.pi, 0, "o", ms=9, mfc="white" if b == 0 else PALETTE[2], mec=PALETTE[2], mew=2, zorder=5)
        for k in (-1, 1):
            ax.plot(k * np.pi, 0, "X", ms=10, color=PALETTE[1], mec="white", zorder=5)
        lab0 = "centre" if b == 0 else "stable spiral"
        ax.annotate(lab0, (0, 0), (0.35, -0.75), fontsize=9.5, color=INK,
                    bbox=dict(fc="white", ec="none", alpha=0.85))
        ax.annotate("saddle", (np.pi, 0), (np.pi + 0.25, 0.35), fontsize=9.5, color=INK,
                    bbox=dict(fc="white", ec="none", alpha=0.85))
        ax.set_xlim(th[0], th[-1]); ax.set_ylim(w[0], w[-1])
        ax.set_xticks(np.arange(-2, 3) * np.pi)
        ax.set_xticklabels(["−2π", "−π", "0", "π", "2π"])
        ax.set_xlabel("angle θ (rad)"); ax.set_ylabel("angular velocity ω = θ'")
        ax.set_title(title)
    if True:
        # one damped trajectory that swings over the top once and then gets captured
        sol = solve_ivp(lambda t, u: [u[1], -np.sin(u[0]) - 0.25 * u[1]], (0, 40), [-2 * np.pi + 0.3, 2.9],
                        max_step=0.02)
        axes[1].plot(sol.y[0], sol.y[1], color=PALETTE[4], lw=2.2, label="one trajectory")
        axes[1].legend(loc="upper right", fontsize=9, frameon=True, framealpha=0.9)
    save(fig, "ch83-pendulum-portrait")


def predator_prey():
    a, b, g, d = 1.0, 0.5, 0.75, 0.25
    f = lambda t, u: [a * u[0] - b * u[0] * u[1], d * u[0] * u[1] - g * u[1]]
    fig, (ax, ax2) = new_figure(1, 2, figsize=(11, 4.3))
    X, Y = np.meshgrid(np.linspace(0.05, 11, 260), np.linspace(0.05, 7, 260))
    H = d * X - g * np.log(X) + b * Y - a * np.log(Y)
    U, V = a * X - b * X * Y, d * X * Y - g * Y
    ax.streamplot(X, Y, U, V, color="#c7d2fe", density=1.1, linewidth=0.7, arrowsize=0.8)
    Hs = d * 3 - g * np.log(3) + b * 2 - a * np.log(2)
    ax.contour(X, Y, H, levels=Hs + np.array([0.05, 0.2, 0.5, 1.0, 1.6]), colors=[PALETTE[0]], linewidths=1.8)
    ax.axvline(g / d, color=MUTED, ls="--", lw=1.2)
    ax.axhline(a / b, color=MUTED, ls=":", lw=1.4)
    ax.text(3.15, 6.6, "ẏ = 0  (x = γ/δ)", fontsize=9, color=MUTED, bbox=dict(fc="white", ec="none", alpha=0.9))
    ax.text(8.6, 2.12, "ẋ = 0  (y = α/β)", fontsize=9, color=MUTED, bbox=dict(fc="white", ec="none", alpha=0.9))
    ax.plot(3, 2, "o", ms=9, mfc="white", mec=PALETTE[2], mew=2, zorder=5)
    ax.plot(0, 0, "X", ms=10, color=PALETTE[1], zorder=5, clip_on=False)
    ax.annotate("coexistence: centre", (3, 2), (4.0, 0.6), fontsize=9.5, arrowprops=dict(arrowstyle="-", color=MUTED))
    ax.annotate("extinction: saddle", (0, 0), (0.5, 0.72), fontsize=9.5, arrowprops=dict(arrowstyle="-", color=MUTED),
                bbox=dict(fc="white", ec="none", alpha=0.9))
    ax.set_xlim(0, 11); ax.set_ylim(0, 7)
    ax.set_xlabel("prey x"); ax.set_ylabel("predators y")
    ax.set_title("closed orbits = level curves of H(x, y)")
    sol = solve_ivp(f, (0, 30), [6, 2], max_step=0.01)
    ax2.plot(sol.t, sol.y[0], color=PALETTE[0], label="prey x(t)")
    ax2.plot(sol.t, sol.y[1], color=PALETTE[1], label="predators y(t)")
    ax2.set_xlabel("time t"); ax2.set_ylabel("population")
    ax2.set_title("the same cycle in time: predators lag prey")
    ax2.legend(loc="upper right", ncol=2)
    ax2.set_ylim(0, 8.5)
    save(fig, "ch83-predator-prey")


def van_der_pol():
    mu = 1.0
    f = lambda t, u: [u[1], mu * (1 - u[0] ** 2) * u[1] - u[0]]
    fig, ax = new_figure(figsize=(7, 5))
    X, Y = np.meshgrid(np.linspace(-4, 4, 200), np.linspace(-4, 4, 200))
    ax.streamplot(X, Y, Y, mu * (1 - X ** 2) * Y - X, color="#d1d5db", density=1.0, linewidth=0.7, arrowsize=0.8)
    s_in = solve_ivp(f, (0, 30), [0.1, 0.0], max_step=0.01)
    s_out = solve_ivp(f, (0, 30), [3.6, 3.2], max_step=0.01)
    s_lc = solve_ivp(f, (60, 80), solve_ivp(f, (0, 60), [2, 0]).y[:, -1], max_step=0.01)
    ax.plot(s_in.y[0], s_in.y[1], color=PALETTE[0], lw=1.4, label="start near the origin (0.1, 0)")
    ax.plot(s_out.y[0], s_out.y[1], color=PALETTE[4], lw=1.4, label="start far out (3.6, 3.2)")
    ax.plot(s_lc.y[0], s_lc.y[1], color=PALETTE[1], lw=3.0, label="limit cycle")
    ax.plot(0, 0, "o", ms=8, mfc="white", mec=PALETTE[1], mew=2, zorder=5)
    ax.annotate("unstable spiral", (0, 0), (0.4, -0.9), fontsize=9.5,
                bbox=dict(fc="white", ec="none", alpha=0.85))
    ax.set_xlim(-4, 4); ax.set_ylim(-4, 4); ax.set_aspect("equal")
    ax.set_xlabel("x"); ax.set_ylabel("x'")
    ax.set_title("Van der Pol, μ = 1: orbits are drawn to one cycle")
    ax.legend(loc="upper left", fontsize=8.5, frameon=True, framealpha=0.9)
    save(fig, "ch83-van-der-pol")


if __name__ == "__main__":
    pendulum_portrait()
    predator_prey()
    van_der_pol()
