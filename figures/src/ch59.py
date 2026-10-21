"""Figures for Chapter 59 — Differential Geometry: Mathematics on Curved Worlds."""
import sys, pathlib; sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
from style import plt, np, PALETTE, INK, MUTED, new_figure, diagram, node, box, arrow, save

BLUE, RED, GREEN, PURPLE, ORANGE = PALETTE[:5]


def sph(lat, lon):
    """Unit-sphere point from latitude/longitude in degrees."""
    la, lo = np.radians(lat), np.radians(lon)
    return np.array([np.cos(la) * np.cos(lo), np.cos(la) * np.sin(lo), np.sin(la)])


def slerp(p, q, n=200):
    w = np.arccos(np.clip(p @ q, -1, 1))
    t = np.linspace(0, 1, n)[:, None]
    return (np.sin((1 - t) * w) * p + np.sin(t * w) * q) / np.sin(w)


# ---------------------------------------------------------------- figure 1: geodesic vs latitude line
fig = plt.figure(figsize=(7.2, 5.6), constrained_layout=True)
ax = fig.add_subplot(projection="3d")
ax.set_box_aspect((1, 1, 1), zoom=1.35)
u = np.linspace(0, 2 * np.pi, 60)
for lat in range(-60, 90, 30):
    la = np.radians(lat)
    ax.plot(np.cos(la) * np.cos(u), np.cos(la) * np.sin(u), np.sin(la) + 0 * u, color="#d1d5db", lw=0.7)
v = np.linspace(-np.pi / 2, np.pi / 2, 40)
for lon in range(0, 360, 30):
    lo = np.radians(lon)
    ax.plot(np.cos(v) * np.cos(lo), np.cos(v) * np.sin(lo), np.sin(v), color="#d1d5db", lw=0.7)

LAT = 40.0
A, B = sph(LAT, -55), sph(LAT, 25)
gc = slerp(A, B)
lons = np.linspace(-55, 25, 200)
lat_line = np.array([sph(LAT, l) for l in lons])
ax.plot(*lat_line.T, color=RED, lw=2.4, label="latitude line (constant heading): not a geodesic")
ax.plot(*gc.T, color=BLUE, lw=2.6, label="great circle: the geodesic, shorter")
ax.scatter(*np.c_[A, B], color=INK, s=30, depthshade=False)
ax.text(*(A * 1.08), "A", color=INK, fontsize=12)
ax.text(*(B * 1.08 + [0, 0, 0.02]), "B", color=INK, fontsize=12)

# tangent plane at A
nA = A
e1 = np.cross([0, 0, 1], nA); e1 /= np.linalg.norm(e1)
e2 = np.cross(nA, e1)
s = np.linspace(-0.45, 0.45, 2)
S, T = np.meshgrid(s, s)
P = A[:, None, None] + e1[:, None, None] * S + e2[:, None, None] * T
ax.plot_surface(P[0], P[1], P[2], color=GREEN, alpha=0.25, edgecolor=GREEN, lw=0.6)
d = (gc[5] - gc[0]); d /= np.linalg.norm(d)
ax.quiver(*A, *(0.4 * d), color=BLUE, lw=2, arrow_length_ratio=0.25)
ax.text(*(A + e2 * 0.5 + e1 * 0.05), r"tangent plane $T_A S^2$", color=GREEN, fontsize=10)
ax.set_xlim(-1, 1); ax.set_ylim(-1, 1); ax.set_zlim(-1, 1)
ax.view_init(elev=22, azim=-15)
ax.set_axis_off()
ax.legend(loc="lower center", fontsize=9.5, bbox_to_anchor=(0.5, -0.02))
ax.set_title("Two routes between points on the same latitude", fontsize=12)
save(fig, "ch59-geodesic-sphere")


# ---------------------------------------------------------------- figure 2: angle sums on three geometries
def poincare_arc(p, q, n=4000):
    """Hyperbolic geodesic between p and q in the Poincare disk, as points."""
    p, q = np.asarray(p, float), np.asarray(q, float)
    if abs(p[0] * q[1] - p[1] * q[0]) < 1e-9:
        return np.linspace(p, q, n)
    pi = p / (p @ p)  # inversion of p in the unit circle lies on the same circle
    M = 2 * np.array([q - p, pi - p])
    rhs = np.array([q @ q - p @ p, pi @ pi - p @ p])
    c = np.linalg.solve(M, rhs)
    r = np.linalg.norm(p - c)
    a0, a1 = np.arctan2(*(p - c)[::-1]), np.arctan2(*(q - c)[::-1])
    da = (a1 - a0 + np.pi) % (2 * np.pi) - np.pi
    t = a0 + np.linspace(0, 1, n) * da
    return np.c_[c[0] + r * np.cos(t), c[1] + r * np.sin(t)]


def angle_at(curve_a, curve_b):
    da = curve_a[1] - curve_a[0]; db = curve_b[1] - curve_b[0]
    return np.degrees(np.arccos(da @ db / np.linalg.norm(da) / np.linalg.norm(db)))


fig = plt.figure(figsize=(11.0, 3.9), constrained_layout=True)

# flat
ax = fig.add_subplot(1, 3, 1)
P = np.array([[0, 0], [4, 0], [1.3, 3.0], [0, 0]])
ax.fill(P[:, 0], P[:, 1], color=BLUE, alpha=0.15); ax.plot(P[:, 0], P[:, 1], color=BLUE)
ax.set_aspect("equal"); ax.axis("off"); ax.set_xlim(-0.5, 4.5); ax.set_ylim(-0.8, 3.6)
ax.set_title("flat plane, $K=0$", fontsize=11)
ax.text(2.0, -0.6, r"$\alpha+\beta+\gamma = 180^\circ$", ha="center", fontsize=11, color=INK)

# sphere octant
ax = fig.add_subplot(1, 3, 2, projection="3d")
ax.set_box_aspect((1, 1, 1))
uu, vv = np.meshgrid(np.linspace(0, 2 * np.pi, 40), np.linspace(0, np.pi, 20))
ax.plot_wireframe(np.cos(uu) * np.sin(vv), np.sin(uu) * np.sin(vv), np.cos(vv), color="#d1d5db", lw=0.5)
N, E1, E2 = np.array([0, 0, 1.0]), np.array([1.0, 0, 0]), np.array([0, 1.0, 0])
for p, q in [(N, E1), (E1, E2), (E2, N)]:
    ax.plot(*slerp(p, q).T, color=RED, lw=2.6)
tt = np.linspace(0, np.pi / 2, 30)
for th in np.linspace(0.05, np.pi / 2 - 0.05, 9):
    ax.plot(np.cos(tt) * np.sin(th), np.sin(tt) * np.sin(th), np.cos(th) + 0 * tt, color=RED, alpha=0.25, lw=1)
ax.view_init(elev=25, azim=35); ax.set_axis_off()
ax.set_xlim(-1, 1); ax.set_ylim(-1, 1); ax.set_zlim(-1, 1)
ax.set_title("sphere, $K>0$", fontsize=11)
ax.text2D(0.5, 0.02, r"$90^\circ+90^\circ+90^\circ = 270^\circ$", transform=ax.transAxes, ha="center", fontsize=11, color=INK)

# hyperbolic (Poincare disk)
ax = fig.add_subplot(1, 3, 3)
t = np.linspace(0, 2 * np.pi, 300)
ax.plot(np.cos(t), np.sin(t), color=MUTED, lw=1.2)
V = [np.array([0.75 * np.cos(a), 0.75 * np.sin(a)]) for a in np.radians([90, 210, 330])]
arcs = [poincare_arc(V[i], V[(i + 1) % 3]) for i in range(3)]
for c in arcs:
    ax.plot(c[:, 0], c[:, 1], color=GREEN, lw=2.4)
poly = np.vstack(arcs)
ax.fill(poly[:, 0], poly[:, 1], color=GREEN, alpha=0.15)
angs = [angle_at(arcs[i], arcs[i - 1][::-1]) for i in range(3)]
ax.set_aspect("equal"); ax.axis("off"); ax.set_xlim(-1.1, 1.1); ax.set_ylim(-1.35, 1.1)
ax.set_title("hyperbolic plane (Poincaré disk), $K<0$", fontsize=11)
ax.text(0, -1.27, rf"$\alpha+\beta+\gamma \approx {sum(angs):.0f}^\circ < 180^\circ$", ha="center", fontsize=11, color=INK)
save(fig, "ch59-angle-sums")
