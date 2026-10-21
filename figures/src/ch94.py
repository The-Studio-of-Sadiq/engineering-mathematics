"""Figures for Chapter 94 — Algebraic Geometry: Solving Equations Through Shapes."""
import sys, pathlib; sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
from style import plt, np, PALETTE, INK, MUTED, new_figure, diagram, node, box, arrow, save

BLUE, RED, GREEN, PURPLE, ORANGE, TEAL, GREY = PALETTE


def curve(ax, f, lim, color=BLUE, lw=2.4, n=900):
    xs = np.linspace(lim[0], lim[1], n)
    ys = np.linspace(lim[2], lim[3], n)
    X, Y = np.meshgrid(xs, ys)
    ax.contour(X, Y, f(X, Y), levels=[0], colors=color, linewidths=lw)


def frame(ax, lim, title):
    ax.set_xlim(lim[0], lim[1]); ax.set_ylim(lim[2], lim[3])
    ax.set_aspect("equal")
    ax.axhline(0, color="#d1d5db", lw=0.8, zorder=0); ax.axvline(0, color="#d1d5db", lw=0.8, zorder=0)
    ax.set_title(title, fontsize=10.5)
    ax.tick_params(labelsize=8)


# ---------------------------------------------------------------- figure 1: a gallery of plane curves
fig, axs = plt.subplots(1, 4, figsize=(10.0, 3.3), constrained_layout=True)
L = (-1.6, 1.6, -1.6, 1.6)
items = [
    (lambda x, y: x**2 + y**2 - 1, r"circle $x^2+y^2=1$", None, "smooth everywhere"),
    (lambda x, y: y**2 - x**3, r"cusp $y^2=x^3$", (0, 0), "singular: one tangent,\ncounted twice"),
    (lambda x, y: y**2 - x**2 * (x + 1), r"node $y^2=x^2(x+1)$", (0, 0), "singular: two\ncrossing branches"),
    (lambda x, y: y**2 - (x**3 - x), r"$y^2=x^3-x$", None, "smooth, two real pieces;\nover $\\mathbb{C}$ a torus"),
]
for ax, (f, title, sing, note) in zip(axs, items):
    curve(ax, f, (-1.8, 1.8, -1.8, 1.8))
    frame(ax, L, title)
    if sing is not None:
        ax.plot(*sing, "o", mfc="white", mec=RED, mew=2.2, ms=9, zorder=5)
    ax.text(0, -1.45, note, ha="center", va="bottom", fontsize=8.5, color=RED if sing else MUTED,
            bbox=dict(fc="white", ec="none", pad=1.5))
# tangent lines of the node
axs[2].plot([-1, 1], [-1, 1], color=ORANGE, lw=1.2, ls="--")
axs[2].plot([-1, 1], [1, -1], color=ORANGE, lw=1.2, ls="--")
axs[2].text(1.02, 0.62, "y = x", color=ORANGE, fontsize=8.5)
axs[2].text(-1.55, 1.12, "y = −x", color=ORANGE, fontsize=8.5)
fig.suptitle("Plane algebraic curves: zero sets of one polynomial in two variables",
             fontsize=11, fontweight="bold")
save(fig, "ch94-curve-gallery")


# ---------------------------------------------------------------- figure 2: Bezout's count
fig, axs = plt.subplots(1, 3, figsize=(10.0, 3.7), constrained_layout=True)

# (a) ellipse and hyperbola: four real points
ax = axs[0]; lim = (-2.6, 2.6, -2.0, 2.0)
curve(ax, lambda x, y: x**2 + 4 * y**2 - 4, (-3, 3, -2.5, 2.5), BLUE)
curve(ax, lambda x, y: x**2 - y**2 - 1, (-3, 3, -2.5, 2.5), GREEN)
xp, yp = np.sqrt(8 / 5), np.sqrt(3 / 5)
for sx in (-1, 1):
    for sy in (-1, 1):
        ax.plot(sx * xp, sy * yp, "o", color=RED, ms=7, zorder=5)
frame(ax, lim, "ellipse ∩ hyperbola: 4 points")
ax.text(0, -1.85, r"$x^2+4y^2=4,\ x^2-y^2=1$", ha="center", fontsize=8.5, color=MUTED)

# (b) two circles: the robot elbow, two real points and two at infinity
ax = axs[1]; lim = (-2.4, 3.4, -2.4, 2.4)
curve(ax, lambda x, y: x**2 + y**2 - 4, (-3, 4, -3, 3), BLUE)
curve(ax, lambda x, y: (x - 2)**2 + (y - 1)**2 - 1, (-3, 4, -3, 3), GREEN)
for (px, py) in [(2.0, 0.0), (1.2, 1.6)]:
    ax.plot(px, py, "o", color=RED, ms=7, zorder=5)
    ax.plot([0, px, 2], [0, py, 1], color=GREY, lw=1.6, zorder=4)
ax.plot(0, 0, "s", color=INK, ms=6, zorder=6); ax.plot(2, 1, "*", color=ORANGE, ms=13, zorder=6)
ax.text(2.15, 1.15, "target", fontsize=8.5, color=ORANGE)
ax.text(-0.15, -0.35, "base", fontsize=8.5, color=INK, ha="right")
frame(ax, lim, "circle ∩ circle: only 2 points?")
ax.text(0.5, -2.25, r"missing two: $(1:\pm j:0)$ at infinity", ha="center", fontsize=8.5,
        color=RED)

# (c) tangency: parabola and circle, a point that counts twice
ax = axs[2]; lim = (-2.0, 2.0, -0.8, 2.4)
curve(ax, lambda x, y: y - x**2, (-2.2, 2.2, -1, 2.6), BLUE)
curve(ax, lambda x, y: x**2 + (y - 1)**2 - 1, (-2.2, 2.2, -1, 2.6), GREEN)
for (px, py) in [(-1, 1), (1, 1)]:
    ax.plot(px, py, "o", color=RED, ms=7, zorder=5)
ax.plot(0, 0, "o", mfc="white", mec=RED, mew=2.4, ms=11, zorder=5)
ax.text(0, -0.55, "tangent: multiplicity 2", ha="center", fontsize=8.5, color=RED)
frame(ax, lim, "parabola ∩ circle: 1 + 1 + 2 = 4")
fig.suptitle("Bézout: two conics meet in 4 points, once you count correctly",
             fontsize=11, fontweight="bold")
save(fig, "ch94-bezout-conics")
