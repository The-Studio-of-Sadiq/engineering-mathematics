"""Chapter 74 figures: gradient-noise fBm heightmap, cellular-automaton caves, an L-system plant."""
import sys, pathlib; sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
from style import plt, np, PALETTE, INK, MUTED, new_figure, diagram, node, box, arrow, save
from matplotlib.colors import LinearSegmentedColormap

BLUE, RED, GREEN, PURPLE, ORANGE, TEAL, GREY = PALETTE

# ---------------------------------------------------------------- 1. gradient noise and fBm
rng = np.random.default_rng(74)
PERM = rng.permutation(256)
PERM = np.concatenate([PERM, PERM])
ANG = rng.random(256) * 2 * np.pi
GRAD = np.stack([np.cos(ANG), np.sin(ANG)], axis=-1)


def fade(t):
    return t * t * t * (t * (6 * t - 15) + 10)


def perlin(x, y):
    xi, yi = np.floor(x).astype(int), np.floor(y).astype(int)
    xf, yf = x - xi, y - yi
    xi, yi = xi & 255, yi & 255

    def dot(ix, iy, dx, dy):
        g = GRAD[PERM[PERM[ix] + iy]]
        return g[..., 0] * dx + g[..., 1] * dy

    n00 = dot(xi, yi, xf, yf)
    n10 = dot(xi + 1, yi, xf - 1, yf)
    n01 = dot(xi, yi + 1, xf, yf - 1)
    n11 = dot(xi + 1, yi + 1, xf - 1, yf - 1)
    u, v = fade(xf), fade(yf)
    return (n00 * (1 - u) + n10 * u) * (1 - v) + (n01 * (1 - u) + n11 * u) * v


def fbm(x, y, octaves, gain=0.5, lac=2.0):
    total, amp, freq = 0.0, 1.0, 1.0
    for _ in range(octaves):
        total = total + amp * perlin(x * freq, y * freq)
        amp *= gain; freq *= lac
    return total


n = 400
xs = np.linspace(0, 6, n, endpoint=False)
X, Y = np.meshgrid(xs, xs)
terrain = LinearSegmentedColormap.from_list("terr", [
    (0.00, "#1e3a8a"), (0.42, "#3b82f6"), (0.47, "#fde68a"), (0.52, "#4ade80"),
    (0.68, "#15803d"), (0.82, "#78716c"), (1.00, "#f8fafc")])
fig, axes = new_figure(1, 3, figsize=(11.0, 4.0))
for ax, k in zip(axes, (1, 3, 7)):
    Z = fbm(X, Y, k)
    lim = np.abs(Z).max()
    ax.imshow(Z, cmap=terrain, vmin=-lim, vmax=lim, origin="lower", interpolation="bilinear")
    ax.set_title(f"{k} octave{'s' if k > 1 else ''}  (gain 1/2)", fontsize=10.5)
    ax.set_xticks([]); ax.set_yticks([]); ax.grid(False)
save(fig, "ch74-noise-heightmap")

# ---------------------------------------------------------------- 2. cellular-automaton caves
def ca_step(G):
    P = np.pad(G, 1, constant_values=1)              # outside the map counts as wall
    S = sum(P[1 + dy:P.shape[0] - 1 + dy, 1 + dx:P.shape[1] - 1 + dx]
            for dy in (-1, 0, 1) for dx in (-1, 0, 1))  # walls in the 3x3 block, centre included
    return (S >= 5).astype(int)


rng = np.random.default_rng(7)
G = (rng.random((60, 90)) < 0.45).astype(int)
snaps = [(0, G.copy())]
for step in range(1, 6):
    G = ca_step(G)
    if step in (1, 5):
        snaps.append((step, G.copy()))
fig, axes = new_figure(1, 3, figsize=(11.0, 2.9))
cave = LinearSegmentedColormap.from_list("cave", ["#fef3c7", "#44403c"])
for ax, (s, g) in zip(axes, snaps):
    ax.imshow(g, cmap=cave, interpolation="nearest")
    ax.set_title(f"step {s}:  {100 * g.mean():.0f}% wall", fontsize=10.5)
    ax.set_xticks([]); ax.set_yticks([]); ax.grid(False)
save(fig, "ch74-cave-automaton")

# ---------------------------------------------------------------- 3. L-system plant
rules = {"X": "F+[[X]-X]-F[-FX]+X", "F": "FF"}
s = "X"
for _ in range(5):
    s = "".join(rules.get(ch, ch) for ch in s)
ang = np.radians(25)
pos, heading, stack, segs = np.array([0.0, 0.0]), np.radians(90 - 10), [], []
for ch in s:
    if ch == "F":
        new = pos + np.array([np.cos(heading), np.sin(heading)])
        segs.append((pos, new, len(stack)))
        pos = new
    elif ch == "+":
        heading += ang
    elif ch == "-":
        heading -= ang
    elif ch == "[":
        stack.append((pos.copy(), heading))
    elif ch == "]":
        pos, heading = stack.pop()
from matplotlib.collections import LineCollection
depth = np.array([d for _, _, d in segs])
lines = [(a, b) for a, b, _ in segs]
greens = plt.get_cmap("Greens")
cols = [greens(0.95 - 0.55 * min(d, 6) / 6) for d in depth]
widths = [2.2 * 0.72 ** min(d, 6) for d in depth]
fig, (a0, a1) = new_figure(1, 2, figsize=(10.0, 5.2), gridspec_kw={"width_ratios": [1.2, 1.5]})
# left: the first three rewrites, each scaled to the same height
for row, n_it in enumerate((1, 2, 3)):
    t = "X"
    for _ in range(n_it):
        t = "".join(rules.get(ch, ch) for ch in t)
    p, hd, st, seg = np.array([0.0, 0.0]), np.radians(80), [], []
    for ch in t:
        if ch == "F":
            q = p + np.array([np.cos(hd), np.sin(hd)]); seg.append((p, q)); p = q
        elif ch == "+": hd += ang
        elif ch == "-": hd -= ang
        elif ch == "[": st.append((p.copy(), hd))
        elif ch == "]": p, hd = st.pop()
    Q = np.array([pt for ab in seg for pt in ab])
    sc = 8.0 / (Q[:, 1].max() - Q[:, 1].min())
    cx = Q[:, 0].mean()
    off = np.array([row * 6.0, 0.0])
    a0.add_collection(LineCollection([((a - [cx, 0]) * sc + off, (b - [cx, 0]) * sc + off) for a, b in seg],
                                     colors=GREEN, linewidths=1.6))
    a0.text(off[0], -1.2, f"n = {n_it}\n{t.count('F')} F", ha="center", va="top", fontsize=9.5, color=INK)
a0.set_xlim(-3, 15); a0.set_ylim(-3.2, 9); a0.set_aspect("equal"); a0.axis("off")
a0.set_title("First rewrites", fontsize=11)
a1.add_collection(LineCollection(lines, colors=cols, linewidths=widths, capstyle="round"))
P = np.array([p for ab in lines for p in ab])
a1.set_xlim(P[:, 0].min() - 2, P[:, 0].max() + 2); a1.set_ylim(P[:, 1].min() - 2, P[:, 1].max() + 2)
a1.set_aspect("equal"); a1.axis("off")
a1.set_title(f"n = 5: {s.count('F')} segments from one 18-symbol rule", fontsize=11)
save(fig, "ch74-lsystem-plant")
print("F count n=5:", s.count("F"), "length", len(s))
print("CA wall fractions", [round(g.mean(), 3) for _, g in snaps])
