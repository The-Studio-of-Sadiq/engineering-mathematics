import sys, pathlib; sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
from style import plt, np, PALETTE, INK, MUTED, new_figure, diagram, node, box, arrow, save
from math import comb
from matplotlib.patches import RegularPolygon, Rectangle

# ---------------------------------------------------------------- Pascal's triangle, odd entries highlighted
ROWS = 16
fig, ax = plt.subplots(figsize=(8.6, 5.6), constrained_layout=True)
ax.set_aspect("equal"); ax.axis("off")
h = np.sqrt(3) / 2
for n in range(ROWS):
    for k in range(n + 1):
        x, y = k - n / 2, -n * h
        v = comb(n, k)
        odd = v % 2 == 1
        ax.add_patch(RegularPolygon((x, y), 6, radius=0.56, orientation=0,
                                    fc=PALETTE[0] if odd else "#f3f4f6",
                                    ec="white", lw=1.2))
        fs = 7.5 if v < 100 else (6.0 if v < 1000 else 4.9)
        ax.text(x, y, str(v), ha="center", va="center", fontsize=fs,
                color="white" if odd else MUTED)
for n in [0, 4, 8, 12]:
    ax.text(-ROWS / 2 - 0.4 + (ROWS - n) * 0 - 0.6, -n * h, f"row {n}", ha="right", va="center",
            fontsize=8.5, color=MUTED)
ax.set_xlim(-ROWS / 2 - 2.6, ROWS / 2 + 0.8)
ax.set_ylim(-(ROWS - 1) * h - 0.8, 0.8)
ax.text(ROWS / 2 - 1.0, -1.2, "blue = odd\ngrey = even", ha="right", va="top", fontsize=9.5, color=INK)
save(fig, "ch28-pascal-parity")

# ---------------------------------------------------------------- pigeonhole and Ramsey R(3,3) = 6
fig, axes = plt.subplots(1, 3, figsize=(11.0, 4.0), constrained_layout=True)
for a in axes:
    a.set_aspect("equal"); a.axis("off")

# panel 1: 5 pigeons, 4 holes
ax = axes[0]
ax.set_xlim(-0.2, 8.2); ax.set_ylim(-0.6, 6.2)
counts = [1, 2, 1, 1]
for i, c in enumerate(counts):
    x0 = 0.2 + 2.0 * i
    ax.add_patch(Rectangle((x0, 0.4), 1.6, 2.6, fc="#f3f4f6", ec=INK, lw=1.5))
    for j in range(c):
        col = PALETTE[1] if c > 1 else PALETTE[0]
        ax.add_patch(plt.Circle((x0 + 0.8, 0.95 + 0.95 * j), 0.36, fc=col, ec="white", lw=1.5))
    ax.text(x0 + 0.8, 0.0, f"box {i + 1}", ha="center", fontsize=9, color=MUTED)
ax.set_title("Pigeonhole: 5 items, 4 boxes", fontsize=11)
ax.text(4.0, 3.9, "some box must hold 2", ha="center", fontsize=10, color=PALETTE[1])

# panel 2: K5 coloured with no monochromatic triangle (pentagon red, pentagram blue)
ax = axes[1]
ang = np.pi / 2 + 2 * np.pi * np.arange(5) / 5
P = np.c_[np.cos(ang), np.sin(ang)] * 2.2
for i in range(5):
    for j in range(i + 1, 5):
        d = (j - i) % 5
        col = PALETTE[1] if d in (1, 4) else PALETTE[0]
        ax.plot(*zip(P[i], P[j]), color=col, lw=2.6, ls="-" if col == PALETTE[1] else "--", zorder=1)
for p in P:
    ax.add_patch(plt.Circle(p, 0.22, fc="white", ec=INK, lw=1.6, zorder=3))
ax.set_xlim(-3, 3); ax.set_ylim(-2.9, 3.0)
ax.set_title("$K_5$: no one-colour triangle", fontsize=11)
ax.text(0, -2.75, "red solid = friends, blue dashed = strangers", ha="center", fontsize=8.5, color=MUTED)

# panel 3: K6 forcing argument
ax = axes[2]
v = np.array([0.0, 0.0])
ang = np.linspace(np.pi * 0.95, np.pi * 0.05, 5)
Q = np.c_[np.cos(ang), np.sin(ang)] * 2.4 + np.array([0, -0.9])
v = np.array([0, -0.9])
v = v + np.array([0, -1.4])
for k, q in enumerate(Q):
    col = PALETTE[1] if k in (0, 2, 4) else PALETTE[0]
    ax.plot(*zip(v, q), color=col, lw=2.6, ls="-" if col == PALETTE[1] else "--", zorder=1)
for a_, b_ in [(0, 2), (2, 4), (0, 4)]:
    ax.plot(*zip(Q[a_], Q[b_]), color=MUTED, lw=1.6, ls=":", zorder=1)
for q in Q:
    ax.add_patch(plt.Circle(q, 0.22, fc="white", ec=INK, lw=1.6, zorder=3))
ax.add_patch(plt.Circle(v, 0.26, fc=INK, ec=INK, zorder=3))
ax.text(v[0] + 0.4, v[1] - 0.15, "you", fontsize=9.5, color=INK, va="center")
cen = Q[[0, 2, 4]].mean(axis=0)
ax.text(cen[0] + 0.45, cen[1] + 0.25, "?", ha="center", va="center", fontsize=13, color=INK)
ax.set_xlim(-3, 3); ax.set_ylim(-2.9, 3.0)
ax.set_title("$K_6$: a triangle is forced", fontsize=11)
ax.text(0, -2.75, "5 edges, 2 colours: 3 share a colour", ha="center", fontsize=8.5, color=MUTED)
save(fig, "ch28-pigeonhole-ramsey")
