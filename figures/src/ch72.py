"""Chapter 72 figures: explicit vs symplectic Euler orbits, bounding volumes and a BVH, tunnelling."""
import sys, pathlib; sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
from style import plt, np, PALETTE, INK, MUTED, new_figure, diagram, node, box, arrow, save
from matplotlib.patches import Rectangle, Polygon, Circle

BLUE, RED, GREEN, PURPLE, ORANGE, TEAL, GREY = PALETTE

# ---------------------------------------------------------------- 1. orbits
def accel(p):
    r = np.linalg.norm(p)
    return -p / r ** 3


def run(method, h, steps, p0=(1.0, 0.0), v0=(0.0, 1.2)):
    p, v = np.array(p0), np.array(v0)
    P, E = [p.copy()], [0.5 * v @ v - 1 / np.linalg.norm(p)]
    for _ in range(steps):
        if method == "explicit":
            p, v = p + h * v, v + h * accel(p)
        else:  # symplectic (semi-implicit) Euler: velocity first, then position with the new velocity
            v = v + h * accel(p)
            p = p + h * v
        P.append(p.copy()); E.append(0.5 * v @ v - 1 / np.linalg.norm(p))
    return np.array(P), np.array(E)


h, steps = 0.02, 3000
Pe, Ee = run("explicit", h, steps)
Ps, Es = run("symplectic", h, steps)
t = np.arange(steps + 1) * h
fig, (a1, a2) = new_figure(1, 2, figsize=(10.0, 4.4), gridspec_kw={"width_ratios": [1, 1.15]})
a1.plot(Pe[:, 0], Pe[:, 1], color=RED, lw=1.1, label="explicit Euler")
a1.plot(Ps[:, 0], Ps[:, 1], color=BLUE, lw=1.4, label="symplectic Euler")
a1.scatter([0], [0], s=120, color=ORANGE, zorder=5, edgecolor=INK)
a1.text(0.08, -0.17, "sun", color=INK, fontsize=9)
a1.scatter([1], [0], s=25, color=INK, zorder=6)
a1.text(1.05, 0.08, "start", color=INK, fontsize=9)
a1.set_aspect("equal"); a1.set_xlim(-5.0, 2.0); a1.set_ylim(-3.0, 3.0)
a1.set_title(f"Same planet, same step h = {h}, t = 0..{steps * h:.0f}")
a1.legend(loc="lower left", fontsize=9)
a2.plot(t, Ee, color=RED, lw=1.4, label="explicit Euler")
a2.plot(t, Es, color=BLUE, lw=1.4, label="symplectic Euler")
a2.axhline(Ee[0], color=MUTED, ls="--", lw=1, label="true energy")
a2.set_xlabel("time t"); a2.set_ylabel("total energy  v²/2 − 1/r")
a2.set_title("Energy: drift versus bounded wobble")
a2.legend(loc="upper left", fontsize=9)
save(fig, "ch72-euler-vs-symplectic")

# ---------------------------------------------------------------- 2. bounding volumes + BVH
fig, (b1, b2) = new_figure(1, 2, figsize=(10.0, 4.4), gridspec_kw={"width_ratios": [1, 1.1]})
shapes = {
    "A": np.array([[1.0, 1.0], [2.6, 1.3], [1.6, 2.6]]),
    "B": np.array([[2.4, 2.0], [3.6, 2.4], [3.1, 3.6], [2.3, 3.1]]),
    "C": np.array([[6.0, 1.0], [7.4, 1.2], [6.9, 2.3]]),
    "D": np.array([[6.6, 3.0], [8.2, 3.3], [7.6, 4.4], [6.8, 4.1]]),
}
cols = {"A": BLUE, "B": GREEN, "C": PURPLE, "D": ORANGE}
for k, P in shapes.items():
    b1.add_patch(Polygon(P, closed=True, fc=cols[k], ec=INK, alpha=0.55, lw=1.2, zorder=3))
    lo, hi = P.min(0), P.max(0)
    b1.add_patch(Rectangle(lo, *(hi - lo), fill=False, ec=cols[k], lw=1.4, ls="--", zorder=4))
    c = P.mean(0); b1.text(*c, k, ha="center", va="center", fontsize=11, weight="bold", color=INK, zorder=5)
for grp, keys, col in (("L", "AB", RED), ("R", "CD", TEAL)):
    pts = np.vstack([shapes[k] for k in keys])
    lo, hi = pts.min(0) - 0.15, pts.max(0) + 0.15
    b1.add_patch(Rectangle(lo, *(hi - lo), fill=False, ec=col, lw=2.0, zorder=2))
    b1.text(lo[0], hi[1] + 0.1, f"box {grp}", color=col, fontsize=9.5)
b1.annotate("A and B: boxes overlap,\nshapes do not\n(narrow phase decides)", xy=(2.45, 2.3), xytext=(2.0, 3.9),
            fontsize=8.8, color=INK, arrowprops=dict(arrowstyle="-|>", color=MUTED, lw=1.2))
b1.set_xlim(0.4, 8.8); b1.set_ylim(0.5, 5.4); b1.set_aspect("equal")
b1.set_xticks([]); b1.set_yticks([]); b1.grid(False)
for sp in b1.spines.values():
    sp.set_visible(False)
b1.set_title("Axis-aligned bounding boxes, nested")

b2.set_xlim(0, 10); b2.set_ylim(0, 6); b2.axis("off")
box(b2, 5, 5.2, "root box", w=2.2, color=INK)
box(b2, 2.6, 3.4, "box L", w=1.8, color=RED)
box(b2, 7.4, 3.4, "box R", w=1.8, color=TEAL)
for x, k in zip((1.4, 3.8, 6.2, 8.6), "ABCD"):
    node(b2, x, 1.6, k, color=cols[k])
arrow(b2, 5, 5.2, 2.6, 3.4, shrink=18); arrow(b2, 5, 5.2, 7.4, 3.4, shrink=18)
for x0, xs in ((2.6, (1.4, 3.8)), (7.4, (6.2, 8.6))):
    for x in xs:
        arrow(b2, x0, 3.4, x, 1.6, shrink=16)
b2.text(5, 0.55, "query 'what touches A?': test root, then L, then B;\nbox R fails, so C and D are never examined",
        ha="center", fontsize=9, color=MUTED)
b2.set_title("Bounding-volume hierarchy", fontsize=12, weight="bold")
save(fig, "ch72-bounding-volumes")

# ---------------------------------------------------------------- 3. tunnelling
fig, ax = new_figure(figsize=(8.5, 3.0))
ax.add_patch(Rectangle((4.85, -1.2), 0.3, 2.4, fc=GREY, ec=INK, alpha=0.6))
ax.text(5.0, 1.38, "thin wall", ha="center", fontsize=9, color=INK)
xs = [1.0, 3.6, 6.2, 8.8]
for i, x in enumerate(xs):
    ax.add_patch(Circle((x, 0), 0.3, fc=BLUE if i < 2 else RED, ec=INK, alpha=0.8, zorder=3))
    ax.text(x, -0.62, f"frame {i}", ha="center", fontsize=9, color=INK)
ax.annotate("", xy=(6.2, 0.0), xytext=(3.6, 0.0), arrowprops=dict(arrowstyle="-|>", color=ORANGE, lw=1.8, ls="--"))
ax.add_patch(Rectangle((3.6, -0.3), 2.6, 0.6, fc=ORANGE, alpha=0.15, ec=ORANGE, ls=":", zorder=1))
ax.text(3.2, 0.95, "swept volume (frame 1 to 2)\nintersects the wall", ha="center", fontsize=8.8, color=ORANGE)
ax.text(1.0, 0.55, "step = v·h", fontsize=9, color=MUTED, ha="center")
ax.set_xlim(0.3, 9.6); ax.set_ylim(-1.3, 1.7); ax.set_aspect("equal")
ax.set_xticks([]); ax.set_yticks([]); ax.grid(False)
for s in ax.spines.values():
    s.set_visible(False)
ax.set_title("Tunnelling: no frame ever sees an overlap, so no collision is reported")
save(fig, "ch72-tunnelling")
