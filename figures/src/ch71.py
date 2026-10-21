import sys, pathlib; sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
from style import plt, np, PALETTE, INK, MUTED, new_figure, diagram, node, box, arrow, save
from matplotlib.patches import Circle, Rectangle

# ---------------------------------------------------------------- figure 1: a tiny ray tracer
W, H = 200, 120
EPS = 1e-4
spheres = [  # centre, radius, colour, reflectivity
    (np.array([-1.1, 0.7, 4.2]), 0.7, np.array([0.86, 0.15, 0.15]), 0.05),
    (np.array([0.55, 0.9, 5.2]), 0.9, np.array([0.85, 0.88, 0.92]), 0.75),   # the mirror ball
    (np.array([1.9, 0.45, 3.6]), 0.45, np.array([0.15, 0.39, 0.92]), 0.05),
]
LIGHT = np.array([-3.0, 5.0, 1.0])
EYE = np.array([0.0, 1.0, -1.0])
SKY_TOP, SKY_LOW = np.array([0.45, 0.62, 0.92]), np.array([0.92, 0.95, 1.0])


def normalize(v):
    return v / np.linalg.norm(v, axis=-1, keepdims=True)


def hit_sphere(o, d, c, r):
    """Smallest t > EPS with |o + t d - c| = r, else inf (d is unit, so a = 1)."""
    oc = o - c
    b = np.sum(oc * d, axis=-1)
    disc = b * b - (np.sum(oc * oc, axis=-1) - r * r)
    sq = np.sqrt(np.maximum(disc, 0))
    t0, t1 = -b - sq, -b + sq
    t = np.where(t0 > EPS, t0, np.where(t1 > EPS, t1, np.inf))
    return np.where(disc >= 0, t, np.inf)


def hit_ground(o, d):
    t = np.where(np.abs(d[..., 1]) > 1e-9, -o[..., 1] / np.where(d[..., 1] == 0, 1, d[..., 1]), np.inf)
    return np.where(t > EPS, t, np.inf)


def nearest(o, d):
    ts = [hit_sphere(o, d, c, r) for c, r, _, _ in spheres] + [hit_ground(o, d)]
    ts = np.stack(ts, axis=-1)
    return ts.min(axis=-1), ts.argmin(axis=-1)


def trace(o, d, depth):
    t, idx = nearest(o, d)
    sky = SKY_LOW + (SKY_TOP - SKY_LOW) * np.clip(d[..., 1:2], 0, 1)
    out = sky.copy()
    hit = np.isfinite(t)
    if not hit.any():
        return out
    p = o + np.where(hit, t, 0.0)[..., None] * d
    n = np.zeros_like(p); col = np.zeros_like(p); refl = np.zeros(t.shape)
    for k, (c, r, cc, rf) in enumerate(spheres):
        m = hit & (idx == k)
        n[m] = (p[m] - c) / r; col[m] = cc; refl[m] = rf
    g = hit & (idx == len(spheres))
    n[g] = [0, 1, 0]
    checker = (np.floor(p[..., 0]) + np.floor(p[..., 2])) % 2
    col[g] = np.where(checker[g, None] > 0, [0.92, 0.92, 0.88], [0.35, 0.35, 0.38])
    refl[g] = 0.08
    L = LIGHT - p
    dist = np.linalg.norm(L, axis=-1)
    L = L / dist[..., None]
    ts, _ = nearest(p + EPS * 10 * n, L)          # shadow ray
    lit = (ts > dist).astype(float)
    lam = np.clip(np.sum(n * L, axis=-1), 0, None) * lit
    rdir = d - 2 * np.sum(d * n, axis=-1, keepdims=True) * n
    spec = np.clip(np.sum(rdir * L, axis=-1), 0, None) ** 60 * lit
    local = col * (0.12 + 0.88 * lam[..., None]) + 0.6 * spec[..., None]
    shade = local
    if depth > 0:
        bounce = trace(p + EPS * 10 * n, normalize(rdir), depth - 1)
        shade = (1 - refl[..., None]) * local + refl[..., None] * bounce
    out[hit] = shade[hit]
    return out


aspect = W / H
xs = np.linspace(-aspect, aspect, W) * 0.62
ys = np.linspace(1, -1, H) * 0.62
X, Y = np.meshgrid(xs, ys)
dirs = normalize(np.stack([X, Y - 0.12, np.ones_like(X)], axis=-1))
origins = np.broadcast_to(EYE, dirs.shape).copy()
img = np.clip(trace(origins, dirs, 2), 0, 1)

fig, ax = new_figure(figsize=(6.6, 4.2))
ax.imshow(img, interpolation="nearest")
ax.set_xticks([0, 50, 100, 150, 199]); ax.set_yticks([0, 40, 80, 119])
ax.grid(False)
ax.set_xlabel("pixel column"); ax.set_ylabel("pixel row")
ax.set_title(f"{W}×{H} pixels, one ray each: shadows, highlights, a mirror ball", fontsize=10.5)
save(fig, "ch71-render")

# ---------------------------------------------------------------- figure 2: ray-sphere intersection
fig, ax = diagram(figsize=(7.4, 4.2))
ax.set_xlim(-0.3, 10.3); ax.set_ylim(-0.2, 6.0)
c, r = np.array([6.6, 2.8]), 1.9
ax.add_patch(Circle(c, r, fc="#eef2ff", ec=PALETTE[0], lw=2))
ax.plot(*c, "o", color=PALETTE[0], ms=4); ax.text(c[0] + 0.12, c[1] - 0.1, "$\\mathbf{c}$", color=PALETTE[0], fontsize=12, va="top")
o = np.array([0.6, 1.2]); d = np.array([1.0, 0.33]); d = d / np.linalg.norm(d)
oc = o - c; b = oc @ d; disc = b * b - (oc @ oc - r * r)
t0, t1 = -b - np.sqrt(disc), -b + np.sqrt(disc)
p0, p1, pm = o + t0 * d, o + t1 * d, o - b * d
ax.plot(*np.c_[o, o + 10 * d], color=INK, lw=1.4)
ax.annotate("", xy=o + 1.3 * d, xytext=o, arrowprops=dict(arrowstyle="-|>", color=PALETTE[1], lw=2.4))
ax.plot(*o, "o", color=INK, ms=6)
ax.text(o[0] - 0.1, o[1] - 0.25, "eye $\\mathbf{o}$", fontsize=11, ha="center", va="top")
ax.text(*(o + 0.6 * d + [0, 0.18]), "$\\hat{\\mathbf{d}}$", color=PALETTE[1], fontsize=12, ha="center", va="bottom")
for p, lab, dy in [(p0, "$t_0$ (enter)", -0.32), (p1, "$t_1$ (exit)", -0.32)]:
    ax.plot(*p, "o", color=PALETTE[4], ms=7, zorder=5)
    ax.text(p[0], p[1] + dy, lab, color=PALETTE[4], fontsize=10, ha="center", va="top")
nrm = (p0 - c) / r
ax.annotate("", xy=p0 + 1.1 * nrm, xytext=p0, arrowprops=dict(arrowstyle="-|>", color=PALETTE[2], lw=2))
ax.text(*(p0 + 1.15 * nrm + [-0.05, 0.05]), "normal $\\hat{\\mathbf{n}}$", color=PALETTE[2], fontsize=10, ha="right")
ax.plot(*np.c_[c, pm], color=MUTED, lw=1.2, ls="--")
ax.plot(*pm, "s", color=MUTED, ms=4)
ax.text(pm[0] - 0.35, pm[1] + 0.5, "closest\napproach", color=MUTED, fontsize=8.5, va="bottom", ha="center")
ax.text(*(0.5 * (c + pm) + [0.12, 0.0]), "$h$", color=MUTED, fontsize=11)
# a ray that misses
o2 = np.array([0.6, 4.3]); d2 = np.array([1.0, 0.12]); d2 = d2 / np.linalg.norm(d2)
ax.plot(*np.c_[o2, o2 + 9.4 * d2], color=MUTED, lw=1.2, ls=":")
ax.text(*(o2 + 9.4 * d2 + [-0.1, 0.15]), "miss: $h > r$, discriminant $< 0$", color=MUTED, fontsize=9, ha="right", va="bottom")
ax.text(0.2, 5.7, "$\\mathbf{p}(t)=\\mathbf{o}+t\\hat{\\mathbf{d}}$,   $|\\mathbf{p}(t)-\\mathbf{c}|=r$"
        "   $\\Rightarrow$   $t^2+2bt+k=0$", fontsize=11, color=INK)
save(fig, "ch71-ray-sphere")

# ---------------------------------------------------------------- figure 3: reflection, refraction, total internal reflection
fig, axes = plt.subplots(1, 2, figsize=(9.6, 4.3), constrained_layout=True)
for ax in axes:
    ax.set_xlim(-2.2, 2.2); ax.set_ylim(-2.0, 2.0); ax.set_aspect("equal"); ax.axis("off")
    ax.add_patch(Rectangle((-2.2, -2.0), 4.4, 2.0, fc="#e0f2fe", ec="none"))
    ax.plot([-2.2, 2.2], [0, 0], color=PALETTE[5], lw=1.6)
    ax.plot([0, 0], [-1.9, 1.9], color=MUTED, lw=1, ls="--")


def ray(ax, ang_from_normal, side, length, color, label, lw=2.2, incoming=False, ls="-"):
    """Draw a ray from the origin making the given angle with the normal; side=+1 above, -1 below."""
    sgn = -1 if incoming else 1
    dx, dy = sgn * np.sin(ang_from_normal), side * np.cos(ang_from_normal)
    end = np.array([dx, dy]) * length
    if incoming:
        ax.annotate("", xy=(0, 0), xytext=end, arrowprops=dict(arrowstyle="-|>", color=color, lw=lw, ls=ls))
    else:
        ax.annotate("", xy=end, xytext=(0, 0), arrowprops=dict(arrowstyle="-|>", color=color, lw=lw, ls=ls))
    ax.text(*(end * 1.08), label, color=color, fontsize=10, ha="center", va="center")


def angle_arc(ax, a0, a1, rad, color, label):
    t = np.linspace(a0, a1, 40)
    ax.plot(rad * np.cos(t), rad * np.sin(t), color=color, lw=1.2)
    tm = (a0 + a1) / 2
    ax.text(1.35 * rad * np.cos(tm), 1.35 * rad * np.sin(tm), label, color=color, fontsize=10.5, ha="center", va="center")


ax = axes[0]
th1 = np.radians(50); th2 = np.arcsin(np.sin(th1) / 1.5)
ray(ax, th1, +1, 1.8, PALETTE[1], "incident", incoming=True)
ray(ax, th1, +1, 1.6, PALETTE[4], "reflected", lw=1.6)
ray(ax, th2, -1, 1.8, PALETTE[0], "refracted")
angle_arc(ax, np.pi / 2, np.pi / 2 + th1, 0.55, PALETTE[1], "$\\theta_1$")
angle_arc(ax, -np.pi / 2, -np.pi / 2 + th2, 0.65, PALETTE[0], "$\\theta_2$")
ax.text(-2.1, 1.85, "air, $n_1 = 1$", fontsize=10, va="top"); ax.text(-2.1, -1.85, "glass, $n_2 = 1.5$", fontsize=10, va="bottom")
ax.text(0.05, 1.95, "normal", fontsize=8.5, color=MUTED, va="top")
ax.set_title(f"Snell: $n_1\\sin\\theta_1 = n_2\\sin\\theta_2$  (50° in, {np.degrees(th2):.1f}° out)", fontsize=10.5)

ax = axes[1]
crit = np.arcsin(1 / 1.5)
# light inside the glass going up: below is glass (n=1.5), above air
ax.text(-2.1, 1.85, "air, $n = 1$", fontsize=10, va="top"); ax.text(-2.1, -1.85, "glass, $n = 1.5$", fontsize=10, va="bottom")
for a, col, lab in [(np.radians(25), PALETTE[2], "25°"), (crit, PALETTE[4], f"{np.degrees(crit):.1f}° (critical)"),
                    (np.radians(60), PALETTE[1], "60°")]:
    end = np.array([-np.sin(a), -np.cos(a)]) * 1.85
    ax.annotate("", xy=(0, 0), xytext=end, arrowprops=dict(arrowstyle="-|>", color=col, lw=1.8))
    ax.text(end[0] - 0.05, end[1] + 0.02, lab, color=col, fontsize=9, ha="right", va="bottom")
    s = 1.5 * np.sin(a)
    if s < 1 - 1e-9:
        out = np.arcsin(s)
        e2 = np.array([np.sin(out), np.cos(out)]) * 1.7
        ax.annotate("", xy=e2, xytext=(0, 0), arrowprops=dict(arrowstyle="-|>", color=col, lw=1.8))
    elif s <= 1 + 1e-9:
        ax.annotate("", xy=(1.9, 0.03), xytext=(0, 0.03), arrowprops=dict(arrowstyle="-|>", color=col, lw=1.8))
        ax.text(1.0, 0.1, "grazes the surface", color=col, fontsize=8.5, ha="center", va="bottom")
    else:
        e2 = np.array([np.sin(a), -np.cos(a)]) * 1.7
        ax.annotate("", xy=e2, xytext=(0, 0), arrowprops=dict(arrowstyle="-|>", color=col, lw=1.8))
        ax.text(e2[0] + 0.05, e2[1], "totally\nreflected", color=col, fontsize=9, ha="left", va="center")
ax.set_title("Glass to air: beyond 41.8° no light escapes", fontsize=10.5)
save(fig, "ch71-refraction")
