import sys, pathlib; sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
from style import plt, np, PALETTE, INK, MUTED, new_figure, diagram, node, box, arrow, save
from itertools import combinations

# ---- figure 1: mug to donut, by interpolating signed-distance functions (one hole throughout)
X, Y = np.meshgrid(np.linspace(-2.2, 2.6, 500), np.linspace(-2.1, 2.0, 420))


def sd_box(x, y, cx, cy, hw, hh, r=0.15):
    qx, qy = np.abs(x - cx) - hw + r, np.abs(y - cy) - hh + r
    return np.hypot(np.maximum(qx, 0), np.maximum(qy, 0)) + np.minimum(np.maximum(qx, qy), 0) - r


def sd_ring(x, y, cx, cy, R, w):
    return np.abs(np.hypot(x - cx, y - cy) - R) - w


mug = np.minimum(sd_box(X, Y, -0.4, 0.0, 1.05, 1.35), sd_ring(X, Y, 0.95, 0.05, 0.62, 0.17))
donut = sd_ring(X, Y, 0.2, 0.0, 1.05, 0.5)
fig, axes = plt.subplots(1, 5, figsize=(11.0, 2.9), constrained_layout=True)
for ax, t in zip(axes, np.linspace(0, 1, 5)):
    f = (1 - t) * mug + t * donut
    ax.contourf(X, Y, f, levels=[-10, 0], colors=["#bfdbfe"])
    ax.contour(X, Y, f, levels=[0], colors=[PALETTE[0]], linewidths=2)
    ax.set_aspect("equal"); ax.axis("off")
    ax.set_title(f"t = {t:.2f}", fontsize=10, color=INK)
axes[0].text(-0.4, -1.95, "mug (solid body + handle)", ha="center", fontsize=9, color=MUTED)
axes[-1].text(0.2, -1.95, "donut (annulus)", ha="center", fontsize=9, color=MUTED)
fig.suptitle("A continuous deformation that never tears or glues: the handle's hole survives every stage",
             fontsize=11, color=INK, weight="bold")
save(fig, "ch58-mug-donut")

# ---- figure 2: loops on a torus
R0, r0 = 2.0, 0.8
fig = plt.figure(figsize=(9.6, 4.2), constrained_layout=True)
ax = fig.add_subplot(projection="3d")
u, v = np.meshgrid(np.linspace(0, 2 * np.pi, 80), np.linspace(0, 2 * np.pi, 40))
T = lambda u, v: ((R0 + r0 * np.cos(v)) * np.cos(u), (R0 + r0 * np.cos(v)) * np.sin(u), r0 * np.sin(v))
ax.plot_surface(*T(u, v), color="#e5e7eb", alpha=0.35, linewidth=0, rstride=1, cstride=1, shade=True)
s = np.linspace(0, 2 * np.pi, 300)
ax.plot(*T(np.full_like(s, -0.6), s), color=PALETTE[1], lw=2.6, label="a: around the tube (meridian)")
ax.plot(*T(s, np.full_like(s, np.pi / 2)), color=PALETTE[0], lw=2.6, label="b: around the hole (longitude)")
uc, vc = 1.9 + 0.35 * np.cos(s), np.pi / 2 + 0.45 * np.sin(s)
ax.plot(*T(uc, vc), color=PALETTE[2], lw=2.6, label="c: contractible (shrinks to a point)")
ax.plot(*T(2 * s, 3 * s + 0.4), color=PALETTE[3], lw=1.4, ls="--", label="2b + 3a: winds twice and three times")
ax.set_box_aspect((1, 1, 0.38)); ax.set_axis_off()
ax.view_init(elev=38, azim=-60)
ax.legend(loc="center left", fontsize=9.5, bbox_to_anchor=(0.92, 0.5))
save(fig, "ch58-torus-loops")

# ---- figure 3: persistence barcode of a noisy circle (Vietoris-Rips, GF(2) reduction)
rng = np.random.default_rng(7)
n = 40
th = rng.uniform(0, 2 * np.pi, n)
P = np.c_[np.cos(th), np.sin(th)] + rng.normal(0, 0.08, (n, 2))
D = np.linalg.norm(P[:, None] - P[None], axis=2)
cmax = 2.2
edges = sorted([(D[i, j], (i, j)) for i, j in combinations(range(n), 2) if D[i, j] <= cmax])
tris = sorted([(max(D[i, j], D[i, k], D[j, k]), (i, j, k)) for i, j, k in combinations(range(n), 3)
               if max(D[i, j], D[i, k], D[j, k]) <= cmax])
simplices = [(0.0, (i,)) for i in range(n)] + edges + tris
simplices.sort(key=lambda s: (s[0], len(s[1])))
idx = {s[1]: k for k, s in enumerate(simplices)}
cols, low_of, bars = [], {}, {0: [], 1: []}
for k, (val, sim) in enumerate(simplices):
    col = set(idx[f] for f in combinations(sim, len(sim) - 1)) if len(sim) > 1 else set()
    while col and max(col) in low_of:
        col ^= cols[low_of[max(col)]]
    cols.append(col)
    if col:
        low_of[max(col)] = k
        birth = simplices[max(col)][0]
        if val - birth > 1e-9:
            bars[len(sim) - 2].append((birth, val))
paired = set(low_of) | set(low_of.values())
for k, (val, sim) in enumerate(simplices):
    if k not in paired and len(sim) <= 2:
        bars[len(sim) - 1].append((val, cmax))

fig, (a1, a2) = plt.subplots(1, 2, figsize=(10.0, 4.4), constrained_layout=True,
                             gridspec_kw={"width_ratios": [1, 1.5]})
eps = 0.8
for p in P:
    a1.add_patch(plt.Circle(p, eps / 2, fc=PALETTE[0], alpha=0.12, ec="none"))
for d, (i, j) in edges:
    if d <= eps:
        a1.plot(*P[[i, j]].T, color=MUTED, lw=0.8)
a1.plot(*P.T, "o", color=PALETTE[0], ms=5)
a1.set_aspect("equal"); a1.set_xlim(-1.5, 1.5); a1.set_ylim(-1.5, 1.5)
a1.set_title(f"40 noisy circle samples, ε = {eps}", fontsize=10.5)
rows = sorted(bars[0], key=lambda b: -(b[1] - b[0]))
h1 = sorted(bars[1], key=lambda b: -(b[1] - b[0]))
yk = 0
for b, d in h1:
    a2.plot([b, d], [yk, yk], color=PALETTE[1], lw=3 if d - b > 0.3 else 1.6, solid_capstyle="butt"); yk += 1
a2.text(1.68, 0.0, "H1: one long bar\n= one loop", color=PALETTE[1], fontsize=9.5, va="center")
yk += 1
y0 = yk
for b, d in rows:
    a2.plot([b, min(d, cmax)], [yk, yk], color=PALETTE[0], lw=1.6, solid_capstyle="butt"); yk += 1
a2.text(1.0, y0 + 8, "H0: one bar survives\n= one connected piece", color=PALETTE[0], fontsize=9.5, va="center")
a2.axvline(eps, color=MUTED, ls=":", lw=1.2)
a2.text(eps + 0.02, yk - 1, "ε", color=MUTED, fontsize=10)
a2.set_yticks([]); a2.set_xlim(0, cmax); a2.invert_yaxis()
a2.set_xlabel("scale ε (ball diameter)")
a2.set_title("Persistence barcode (Vietoris–Rips filtration)", fontsize=10.5)
save(fig, "ch58-persistence-barcode")
print("H1 bars:", [(round(b, 3), round(d, 3)) for b, d in h1])
print("H0 longest finite:", [round(d - b, 3) for b, d in rows[:4]])
