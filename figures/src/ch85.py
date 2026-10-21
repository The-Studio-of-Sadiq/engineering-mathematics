"""Chapter 85 figures: logistic-map bifurcation diagram, normal-form bifurcation diagrams, subcritical hysteresis."""
import sys, pathlib; sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
from style import plt, np, PALETTE, INK, MUTED, new_figure, diagram, node, box, arrow, save

STABLE = dict(color=PALETTE[0], lw=2.4)
UNSTABLE = dict(color=PALETTE[1], lw=2.0, ls="--")


def logistic_diagram():
    r = np.linspace(2.8, 4.0, 2400)
    x = np.full_like(r, 0.5)
    for _ in range(1000):
        x = r * x * (1 - x)
    fig, ax = new_figure(figsize=(10, 5.2))
    for _ in range(300):
        x = r * x * (1 - x)
        ax.plot(r, x, ",", color=INK, alpha=0.25, rasterized=True)
    marks = [(3.0, "r₁ = 3"), (3.449490, "r₂ = 3.4495"), (3.544090, "r₃"), (3.569946, "r∞ ≈ 3.5699")]
    for i, (rv, lab) in enumerate(marks):
        ax.axvline(rv, color=PALETTE[1] if i < 3 else PALETTE[3], lw=1.0, ls=":")
        ax.text(rv + 0.006, 0.03 + 0.06 * (i % 2), lab, fontsize=9, color=PALETTE[1] if i < 3 else PALETTE[3],
                bbox=dict(fc="white", ec="none", alpha=0.85))
    rr = np.linspace(2.8, 3.0, 50)
    ax.plot(rr, 1 - 1 / rr, color=PALETTE[0], lw=2)
    ax.annotate("one stable fixed point\nx* = 1 − 1/r", (2.88, 1 - 1 / 2.88), (2.82, 0.88), fontsize=9,
                arrowprops=dict(arrowstyle="-", color=MUTED), bbox=dict(fc="white", ec="none", alpha=0.9))
    ax.annotate("period-3 window\nr = 1 + √8 = 3.8284", (3.835, 0.5), (3.62, 0.25), fontsize=9,
                arrowprops=dict(arrowstyle="->", color=MUTED), bbox=dict(fc="white", ec="none", alpha=0.9))
    ax.set_xlim(2.8, 4.0); ax.set_ylim(0, 1)
    ax.set_xlabel("parameter r"); ax.set_ylabel("long-run values of $x_n$")
    ax.set_title("Logistic map $x_{n+1} = r x_n (1 - x_n)$: period doubling to chaos")
    ax.grid(False)
    save(fig, "ch85-logistic-bifurcation")


def flow_arrows(ax, mus, xs, f, scale=0.18):
    for m in mus:
        for xv in xs:
            v = f(m, xv)
            if abs(v) < 1e-3:
                continue
            ax.annotate("", (m, xv + np.sign(v) * scale), (m, xv), arrowprops=dict(arrowstyle="-|>", color=MUTED, lw=1.0))


def normal_forms():
    fig, axes = new_figure(2, 2, figsize=(10, 7.6))
    mu = np.linspace(-1, 1, 400)
    xs = np.linspace(-1.2, 1.2, 9)
    mus = [-0.75, -0.25, 0.3, 0.75]

    ax = axes[0, 0]
    m = mu[mu >= 0]
    ax.plot(m, np.sqrt(m), **STABLE, label="stable")
    ax.plot(m, -np.sqrt(m), **UNSTABLE, label="unstable")
    flow_arrows(ax, mus, xs, lambda a, x: a - x * x)
    ax.set_title("saddle-node:  ẋ = μ − x²")
    ax.legend(loc="upper left", fontsize=9, frameon=True, framealpha=0.9)

    ax = axes[0, 1]
    ax.plot(mu[mu <= 0], 0 * mu[mu <= 0], **STABLE); ax.plot(mu[mu >= 0], 0 * mu[mu >= 0], **UNSTABLE)
    ax.plot(mu[mu >= 0], mu[mu >= 0], **STABLE); ax.plot(mu[mu <= 0], mu[mu <= 0], **UNSTABLE)
    flow_arrows(ax, mus, xs, lambda a, x: a * x - x * x)
    ax.set_title("transcritical:  ẋ = μx − x²")

    ax = axes[1, 0]
    ax.plot(mu[mu <= 0], 0 * mu[mu <= 0], **STABLE); ax.plot(mu[mu >= 0], 0 * mu[mu >= 0], **UNSTABLE)
    m = mu[mu >= 0]
    ax.plot(m, np.sqrt(m), **STABLE); ax.plot(m, -np.sqrt(m), **STABLE)
    flow_arrows(ax, mus, xs, lambda a, x: a * x - x ** 3)
    ax.set_title("supercritical pitchfork:  ẋ = μx − x³")

    ax = axes[1, 1]
    ax.plot(mu[mu <= 0], 0 * mu[mu <= 0], color=PALETTE[0], lw=2.4, label="stable equilibrium")
    ax.plot(mu[mu >= 0], 0 * mu[mu >= 0], **UNSTABLE)
    ax.plot(m, np.sqrt(m), color=PALETTE[2], lw=2.4, label="stable limit cycle, radius √μ")
    ax.plot(m, -np.sqrt(m), color=PALETTE[2], lw=2.4)
    ax.fill_between(m, -np.sqrt(m), np.sqrt(m), color=PALETTE[2], alpha=0.08)
    ax.set_title("Hopf:  ṙ = μr − r³,  θ̇ = ω")
    ax.set_ylabel("x = r cos θ  (cycle spans ±r)")
    ax.legend(loc="upper left", fontsize=8.5, frameon=True, framealpha=0.9)

    for a in axes.ravel():
        a.set_xlim(-1, 1); a.set_ylim(-1.3, 1.3)
        a.axvline(0, color="#d1d5db", lw=0.8, zorder=0)
        a.set_xlabel("parameter μ")
    for a in (axes[0, 0], axes[0, 1], axes[1, 0]):
        a.set_ylabel("equilibrium x*")
    save(fig, "ch85-normal-forms")


def hysteresis():
    fig, ax = new_figure(figsize=(7.4, 4.6))
    # x' = mu x + x^3 - x^5 : nonzero equilibria x^2 = (1 +- sqrt(1 + 4 mu)) / 2
    mu = np.linspace(-0.25, 0.6, 500)
    up = np.sqrt((1 + np.sqrt(1 + 4 * mu)) / 2)
    mu2 = np.linspace(-0.25, 0.0, 300)
    lo = np.sqrt((1 - np.sqrt(1 + 4 * mu2)) / 2)
    ax.plot(np.linspace(-0.6, 0, 50), 0 * np.linspace(-0.6, 0, 50), **STABLE, label="stable")
    ax.plot(np.linspace(0, 0.6, 50), 0 * np.linspace(0, 0.6, 50), **UNSTABLE, label="unstable")
    ax.plot(mu, up, **STABLE); ax.plot(mu2, lo, **UNSTABLE)
    ax.plot(mu, -up, color=PALETTE[0], lw=2.4, alpha=0.35); ax.plot(mu2, -lo, color=PALETTE[1], lw=2.0, ls="--", alpha=0.35)
    kw = dict(arrowprops=dict(arrowstyle="-|>", color=PALETTE[4], lw=2.2))
    ax.annotate("", (-0.02, 0.02), (-0.45, 0.02), **kw)
    ax.annotate("", (0.0, 0.98), (0.0, 0.06), **kw)
    ax.annotate("", (-0.24, 0.73), (0.35, 1.13), arrowprops=dict(arrowstyle="-|>", color=PALETTE[4], lw=2.2,
                                                                connectionstyle="arc3,rad=0.15"))
    ax.annotate("", (-0.25, 0.04), (-0.25, 0.68), **kw)
    ax.text(-0.44, -0.13, "1. increase μ slowly", fontsize=9, color=PALETTE[4])
    ax.text(0.02, 0.45, "2. jump at μ = 0", fontsize=9, color=PALETTE[4])
    ax.text(0.05, 1.27, "3. now decrease μ: stays up", fontsize=9, color=PALETTE[4])
    ax.text(-0.58, 0.40, "4. falls at\nμ = −1/4", fontsize=9, color=PALETTE[4])
    ax.set_xlim(-0.6, 0.6); ax.set_ylim(-1.35, 1.45)
    ax.set_xlabel("parameter μ"); ax.set_ylabel("equilibrium x*")
    ax.set_title("Subcritical pitchfork ẋ = μx + x³ − x⁵: a hysteresis loop")
    ax.legend(loc="lower left", fontsize=9, frameon=True, framealpha=0.9)
    save(fig, "ch85-hysteresis")


if __name__ == "__main__":
    logistic_diagram()
    normal_forms()
    hysteresis()
