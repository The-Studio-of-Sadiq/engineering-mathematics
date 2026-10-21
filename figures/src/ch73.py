"""Chapter 73 figures: escape-time Mandelbrot set, Koch snowflake iterations, box-counting dimension fit."""
import sys, pathlib; sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
from style import plt, np, PALETTE, INK, MUTED, new_figure, diagram, node, box, arrow, save

BLUE, RED, GREEN, PURPLE, ORANGE, TEAL, GREY = PALETTE

# ---------------------------------------------------------------- 1. Mandelbrot by escape time
W, H, MAXIT = 900, 600, 200
re = np.linspace(-2.25, 0.75, W)
im = np.linspace(-1.0, 1.0, H)
c = re[None, :] + 1j * im[:, None]
z = np.zeros_like(c)
count = np.full(c.shape, MAXIT, dtype=float)
alive = np.ones(c.shape, dtype=bool)
for n in range(MAXIT):
    z[alive] = z[alive] ** 2 + c[alive]
    esc = alive & (np.abs(z) > 2)
    count[esc] = n
    alive &= ~esc
shade = np.where(count == MAXIT, np.nan, np.log1p(count))
fig, ax = new_figure(figsize=(8.0, 5.6))
cmap = plt.get_cmap("magma_r").copy(); cmap.set_bad(INK)
ax.imshow(shade, extent=(re[0], re[-1], im[0], im[-1]), origin="lower", cmap=cmap, interpolation="bilinear")
ax.set_xlabel("Re c"); ax.set_ylabel("Im c"); ax.grid(False)
ax.set_title(f"Mandelbrot set: black = no escape in {MAXIT} iterations")
ax.annotate("main cardioid", xy=(-0.1, 0.2), xytext=(0.05, 0.75), color="white", fontsize=9,
            arrowprops=dict(arrowstyle="-|>", color="white", lw=1.1))
ax.annotate("period-2 disc", xy=(-1.0, 0.0), xytext=(-1.6, 0.62), color="white", fontsize=9,
            arrowprops=dict(arrowstyle="-|>", color="white", lw=1.1))
ax.annotate("the cusp, c = 1/4", xy=(0.25, 0.0), xytext=(0.52, -0.88), color=INK, fontsize=9, ha="center",
            arrowprops=dict(arrowstyle="-|>", color="#facc15", lw=1.4))
save(fig, "ch73-mandelbrot")

# ---------------------------------------------------------------- 2. Koch snowflake iterations
def koch_step(P):
    out = []
    for a, b in zip(P[:-1], P[1:]):
        d = (b - a) / 3
        p1, p3 = a + d, a + 2 * d
        p2 = p1 + d * np.exp(1j * np.pi / 3)       # bump points outward for a clockwise triangle
        out += [a, p1, p2, p3]
    out.append(P[-1])
    return np.array(out)


tri = np.exp(1j * (np.pi / 2 - 2 * np.pi * np.arange(4) / 3))   # clockwise, closed
fig, axes = new_figure(1, 4, figsize=(11.0, 3.2))
P = tri
for k, ax in enumerate(axes):
    ax.fill(P.real, P.imag, color="#e5edff", ec=BLUE, lw=1.2)
    n_seg = 3 * 4 ** k
    ax.set_title(f"n = {k}:  {n_seg} sides", fontsize=10.5)
    ax.text(0, -1.45, f"perimeter x {(4 / 3) ** k:.3f}\narea x {1 + 0.6 * (1 - (4 / 9) ** k):.3f}",
            ha="center", va="top", fontsize=9, color=INK)
    ax.set_aspect("equal"); ax.set_xlim(-1.3, 1.3); ax.set_ylim(-2.1, 1.25); ax.axis("off")
    P = koch_step(P)
save(fig, "ch73-koch-iterations")

# ---------------------------------------------------------------- 3. box counting on the Koch curve
seg = np.array([0 + 0j, 1 + 0j])
for _ in range(8):
    a, b = seg[:-1], seg[1:]
    d = (b - a) / 3
    p1, p3 = a + d, a + 2 * d
    p2 = p1 + d * np.exp(1j * np.pi / 3)
    seg = np.column_stack([a, p1, p2, p3]).ravel()
    seg = np.append(seg, b[-1])
# densify so that every box the curve passes through contains a sample point
t = np.linspace(0, 1, 6)[:-1]
pts = (seg[:-1, None] + (seg[1:] - seg[:-1])[:, None] * t[None, :]).ravel()
eps = 2.0 ** -np.arange(2, 11)
Ncount = []
for e in eps:
    keys = np.unique(np.column_stack([np.floor(pts.real / e), np.floor(pts.imag / e)]), axis=0)
    Ncount.append(len(keys))
Ncount = np.array(Ncount)
slope, icept = np.polyfit(np.log(1 / eps), np.log(Ncount), 1)
fig, ax = new_figure(figsize=(6.6, 4.3))
ax.loglog(1 / eps, Ncount, "o", color=BLUE, ms=7, label="box counts  N(ε)")
xx = np.array([1 / eps[0], 1 / eps[-1]])
ax.loglog(xx, np.exp(icept) * xx ** slope, color=RED, lw=1.6, label=f"fit: slope = {slope:.3f}")
ax.loglog(xx, Ncount[0] * (xx / xx[0]) ** 1.0, color=MUTED, ls="--", lw=1.2, label="slope 1 (a smooth curve)")
ax.loglog(xx, Ncount[0] * (xx / xx[0]) ** 2.0, color=MUTED, ls=":", lw=1.4, label="slope 2 (a filled region)")
ax.set_xlabel("1/ε  (boxes per unit length)"); ax.set_ylabel("N(ε): boxes that touch the curve")
ax.set_title(f"Box counting the Koch curve; exact value log 4 / log 3 = {np.log(4) / np.log(3):.3f}", fontsize=10.5)
ax.legend(loc="upper left", fontsize=9)
save(fig, "ch73-box-counting")
print("box-counting slope", slope, "counts", Ncount)
