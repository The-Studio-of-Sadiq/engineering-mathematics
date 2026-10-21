"""Figures for Chapter 75 — Computational Geometry: Algorithms for Shapes."""
import sys, pathlib; sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
from style import plt, np, PALETTE, INK, MUTED, new_figure, diagram, node, box, arrow, save
from scipy.spatial import Delaunay, Voronoi

BLUE, RED, GREEN, PURPLE, ORANGE, TEAL, GREY = PALETTE


def orient(a, b, c):
    return (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])


def graham(P):
    """Andrew's monotone chain variant of Graham's scan; returns hull CCW and the popped points."""
    pts = sorted(map(tuple, P))
    popped = []

    def half(seq):
        h = []
        for p in seq:
            while len(h) >= 2 and orient(h[-2], h[-1], p) <= 0:
                popped.append((h[-2], h[-1], p))
                h.pop()
            h.append(p)
        return h

    lower, upper = half(pts), half(pts[::-1])
    return np.array(lower[:-1] + upper[:-1]), popped


# ---------------------------------------------------------------- figure 1: convex hull
rng = np.random.default_rng(7)
P = rng.normal(size=(40, 2)) * [1.0, 0.7]
H, popped = graham(P)

fig, axs = new_figure(1, 2, figsize=(10.5, 4.4))
ax = axs[0]
ax.set_aspect("equal"); ax.axis("off")
Hc = np.vstack([H, H[:1]])
ax.fill(Hc[:, 0], Hc[:, 1], color=BLUE, alpha=0.10)
ax.plot(Hc[:, 0], Hc[:, 1], color=BLUE, lw=2)
ax.plot(P[:, 0], P[:, 1], "o", color=GREY, ms=4)
ax.plot(H[:, 0], H[:, 1], "o", color=BLUE, ms=6)
ax.set_title(f"convex hull: {len(H)} of {len(P)} points are corners", fontsize=10.5)

# (b) the turn test, zoomed: a left turn kept, a right turn popped
ax = axs[1]
ax.set_aspect("equal"); ax.axis("off")
a, b, c = np.array([0.0, 0.0]), np.array([2.0, 0.6]), np.array([3.4, 2.0])
d = np.array([3.6, -0.2])
ax.annotate("", xy=b, xytext=a, arrowprops=dict(arrowstyle="-|>", color=INK, lw=1.8))
ax.annotate("", xy=c, xytext=b, arrowprops=dict(arrowstyle="-|>", color=GREEN, lw=1.8))
ax.annotate("", xy=d, xytext=b, arrowprops=dict(arrowstyle="-|>", color=RED, lw=1.8, ls="--"))
for p, s, off in [(a, "a", (-0.25, -0.15)), (b, "b", (0.0, 0.3)), (c, "c", (0.15, 0.15)), (d, "c'", (0.2, -0.2))]:
    ax.plot(*p, "o", color=INK, ms=5)
    ax.text(p[0] + off[0], p[1] + off[1], f"${s}$", fontsize=12, color=INK)
ax.text(3.0, 2.55, "left turn: orient > 0\nkeep b", color=GREEN, fontsize=9.5, ha="center")
ax.text(3.2, -0.95, "right turn: orient < 0\npop b", color=RED, fontsize=9.5, ha="center")
ax.text(0.2, 1.9, r"orient$(a,b,c)=\det\,[\,b-a,\ c-a\,]$", fontsize=10, color=INK)
ax.set_xlim(-0.6, 4.8); ax.set_ylim(-1.3, 3.2)
ax.set_title("the only question the scan ever asks", fontsize=10.5)
save(fig, "ch75-convex-hull")

# ---------------------------------------------------------------- figure 2: Delaunay + Voronoi
rng = np.random.default_rng(3)
S = rng.uniform(0, 10, size=(22, 2))
tri = Delaunay(S)
vor = Voronoi(S)

fig, ax = new_figure(figsize=(7.0, 6.0))
ax.set_aspect("equal"); ax.grid(False)
for sp in ax.spines.values():
    sp.set_visible(False)
ax.set_xticks([]); ax.set_yticks([])
ax.triplot(S[:, 0], S[:, 1], tri.simplices, color=BLUE, lw=1.2)
center = S.mean(axis=0)
for (i, j), rv in zip(vor.ridge_points, vor.ridge_vertices):
    rv = np.asarray(rv)
    if np.all(rv >= 0):
        seg = vor.vertices[rv]
    else:
        v = vor.vertices[rv[rv >= 0][0]]
        t = S[j] - S[i]; t /= np.linalg.norm(t)
        n = np.array([-t[1], t[0]])
        mid = (S[i] + S[j]) / 2
        n = n if np.dot(mid - center, n) > 0 else -n
        seg = np.array([v, v + 20 * n])
    ax.plot(seg[:, 0], seg[:, 1], color=RED, lw=1.4, ls="--")
ax.plot(S[:, 0], S[:, 1], "o", color=INK, ms=5, zorder=4)
# highlight one triangle and its empty circumcircle
k = 10
A, B, C = S[tri.simplices[k]]
D = 2 * (A[0] * (B[1] - C[1]) + B[0] * (C[1] - A[1]) + C[0] * (A[1] - B[1]))
ux = ((A @ A) * (B[1] - C[1]) + (B @ B) * (C[1] - A[1]) + (C @ C) * (A[1] - B[1])) / D
uy = ((A @ A) * (C[0] - B[0]) + (B @ B) * (A[0] - C[0]) + (C @ C) * (B[0] - A[0])) / D
R = np.hypot(A[0] - ux, A[1] - uy)
ax.fill(*np.array([A, B, C]).T, color=BLUE, alpha=0.2)
ax.add_patch(plt.Circle((ux, uy), R, fill=False, ec=GREEN, lw=2))
ax.plot(ux, uy, "s", color=GREEN, ms=6, zorder=5)
ax.annotate("Voronoi vertex =\ncircumcentre", xy=(ux, uy), xytext=(ux - 4.2, uy + 1.6), fontsize=9, color=GREEN,
            bbox=dict(fc="white", ec=GREEN, lw=1), arrowprops=dict(arrowstyle="-|>", color=GREEN), zorder=6)
ax.set_xlim(-0.5, 10.5); ax.set_ylim(-0.5, 10.5)
ax.plot([], [], color=BLUE, label="Delaunay triangulation")
ax.plot([], [], color=RED, ls="--", label="Voronoi diagram")
ax.plot([], [], color=GREEN, label="an empty circumcircle")
ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.01), ncol=3, fontsize=9)
ax.set_title("two drawings of one structure", fontsize=11)
save(fig, "ch75-delaunay-voronoi")

# ---------------------------------------------------------------- figure 3: quadtree
rng = np.random.default_rng(11)
cl1 = rng.normal([2.5, 7.0], 0.5, size=(45, 2))
cl2 = rng.normal([7.0, 3.0], 0.9, size=(30, 2))
sc = rng.uniform(0, 10, size=(12, 2))
Q = np.vstack([cl1, cl2, sc])
Q = Q[(Q[:, 0] > 0) & (Q[:, 0] < 10) & (Q[:, 1] > 0) & (Q[:, 1] < 10)]
CAP = 4
depth_count = {}


def build(ax, pts, x0, y0, s, depth):
    depth_count[depth] = depth_count.get(depth, 0) + 1
    ax.add_patch(plt.Rectangle((x0, y0), s, s, fill=False, ec=GREY, lw=max(1.6 - 0.25 * depth, 0.5)))
    if len(pts) <= CAP or depth >= 7:
        return
    h = s / 2
    for dx in (0, 1):
        for dy in (0, 1):
            m = ((pts[:, 0] >= x0 + dx * h) & (pts[:, 0] < x0 + (dx + 1) * h) &
                 (pts[:, 1] >= y0 + dy * h) & (pts[:, 1] < y0 + (dy + 1) * h))
            build(ax, pts[m], x0 + dx * h, y0 + dy * h, h, depth + 1)


fig, axs = new_figure(1, 2, figsize=(10.5, 5.0), gridspec_kw=dict(width_ratios=[1, 1]))
ax = axs[0]
ax.set_aspect("equal"); ax.axis("off")
build(ax, Q, 0, 0, 10, 0)
ax.plot(Q[:, 0], Q[:, 1], "o", color=BLUE, ms=3, zorder=3)
qc, qr = np.array([6.0, 4.0]), 1.6
ax.add_patch(plt.Circle(qc, qr, fill=False, ec=RED, lw=2, zorder=4))
ax.text(qc[0], qc[1] + qr + 0.25, "range query", color=RED, ha="center", fontsize=9,
        bbox=dict(fc="white", ec="none", pad=1), zorder=5)
ax.set_title(f"quadtree, at most {CAP} points per leaf ({len(Q)} points)", fontsize=10.5)

# (b) cells per depth
ax = axs[1]
ds = sorted(depth_count)
ax.bar(ds, [depth_count[d] for d in ds], color=TEAL)
ax.plot(ds, [4 ** d for d in ds], "o--", color=MUTED, label=r"full tree: $4^d$")
ax.set_yscale("log")
ax.set_xlabel("depth $d$"); ax.set_ylabel("number of cells")
ax.set_title("the tree grows only where points crowd", fontsize=10.5)
ax.legend()
save(fig, "ch75-quadtree")
print("hull", len(H), "of", len(P), "quadtree cells", depth_count)
