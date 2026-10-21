import sys, pathlib; sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
from style import plt, np, PALETTE, INK, MUTED, new_figure, diagram, node, box, arrow, save
from matplotlib.patches import Polygon, FancyArrowPatch
from matplotlib.colors import ListedColormap

# Square corners at angles 45, 135, 225, 315 degrees, labelled 1..4 counter-clockwise.
ANG = np.deg2rad([45, 135, 225, 315])
CORNERS = np.c_[np.cos(ANG), np.sin(ANG)]


def rot(k):
    t = k * np.pi / 2
    return np.array([[np.cos(t), -np.sin(t)], [np.sin(t), np.cos(t)]])


def refl(phi):
    """Reflection across the line through the origin at angle phi."""
    c, s = np.cos(2 * phi), np.sin(2 * phi)
    return np.array([[c, s], [s, -c]])


# The eight symmetries as 2x2 matrices, with names and (for reflections) mirror angle.
ELEMS = [("e", rot(0), None), ("r", rot(1), None), ("r²", rot(2), None), ("r³", rot(3), None),
         ("h", refl(0.0), 0.0), ("v", refl(np.pi / 2), np.pi / 2),
         ("d", refl(np.pi / 4), np.pi / 4), ("d′", refl(-np.pi / 4), -np.pi / 4)]


def perm(M):
    """Where each corner label ends up: returns label at each corner position after applying M."""
    img = CORNERS @ M.T
    out = [0] * 4
    for i, p in enumerate(img):
        j = int(np.argmin(np.linalg.norm(CORNERS - p, axis=1)))
        out[j] = i + 1
    return out


# ---- figure 1: the eight symmetries of a square
fig, axes = plt.subplots(2, 4, figsize=(9.6, 5.4), constrained_layout=True)
cols = [PALETTE[0]] * 4 + [PALETTE[1]] * 4
for ax, (name, M, phi), c in zip(axes.flat, ELEMS, cols):
    ax.set_xlim(-1.6, 1.6); ax.set_ylim(-1.75, 1.6); ax.set_aspect("equal"); ax.axis("off")
    labels = perm(M)
    ax.add_patch(Polygon(CORNERS, closed=True, fc="#eef2ff" if c == PALETTE[0] else "#fef2f2", ec=c, lw=2))
    # a little flag so the orientation is visible
    flag = np.array([[0.0, 0.0], [0.4, 0.0], [0.4, 0.2]]) + [-0.2, 0.18]
    ax.add_patch(Polygon(flag @ M.T, closed=True, fc=c, ec=c, alpha=0.8))
    for (x, y), lab in zip(CORNERS, labels):
        ax.text(1.28 * x, 1.28 * y, str(lab), ha="center", va="center", fontsize=11, color=INK, weight="bold")
    if phi is not None:
        L = 1.35 if abs(np.sin(2 * phi)) < 0.5 else 1.08
        t = np.linspace(-L, L, 2)
        ax.plot(t * np.cos(phi), t * np.sin(phi), ls="--", color=c, lw=1.4)
    elif name != "e":
        k = {"r": 1, "r²": 2, "r³": 3}[name]
        th = np.linspace(-0.2, k * np.pi / 2 - 0.35, 60)
        ax.plot(0.55 * np.cos(th), 0.55 * np.sin(th) - 0.15, color=c, lw=1.4)
        ax.annotate("", xy=(0.55 * np.cos(th[-1] + 0.25), 0.55 * np.sin(th[-1] + 0.25) - 0.15),
                    xytext=(0.55 * np.cos(th[-1]), 0.55 * np.sin(th[-1]) - 0.15),
                    arrowprops=dict(arrowstyle="-|>", color=c, lw=1.4, mutation_scale=11))
    desc = {"e": "identity", "r": "rotate 90°", "r²": "rotate 180°", "r³": "rotate 270°",
            "h": "flip, horizontal axis", "v": "flip, vertical axis", "d": "flip, diagonal", "d′": "flip, anti-diagonal"}[name]
    ax.text(0, -1.62, f"{name}: {desc}", ha="center", va="center", fontsize=9.5, color=c)
fig.suptitle("The dihedral group D4: eight ways to put a square back in its footprint", fontsize=11.5, color=INK, weight="bold")
save(fig, "ch56-square-symmetries")

# ---- figure 2: Cayley table as a heat map
names = [n for n, _, _ in ELEMS]
mats = [M for _, M, _ in ELEMS]


def index_of(M):
    return int(np.argmin([np.abs(M - X).sum() for X in mats]))


T = np.array([[index_of(a @ b) for b in mats] for a in mats])  # row a, column b: a∘b (b first)
cmap = ListedColormap(["#1e3a8a", "#2563eb", "#60a5fa", "#bfdbfe", "#7f1d1d", "#dc2626", "#f87171", "#fecaca"])
fig, ax = new_figure(figsize=(6.4, 5.6))
ax.imshow(T, cmap=cmap, vmin=-0.5, vmax=7.5)
ax.grid(False)
for i in range(8):
    for j in range(8):
        k = T[i, j]
        ax.text(j, i, names[k], ha="center", va="center", fontsize=11,
                color="white" if k in (0, 1, 4, 5) else INK, weight="bold")
ax.set_xticks(range(8)); ax.set_xticklabels(names, fontsize=11)
ax.set_yticks(range(8)); ax.set_yticklabels(names, fontsize=11)
ax.xaxis.tick_top(); ax.xaxis.set_label_position("top")
ax.set_xlabel("b  (applied first)"); ax.set_ylabel("a  (applied second)")
ax.axhline(3.5, color="white", lw=3); ax.axvline(3.5, color="white", lw=3)
for s in ax.spines.values():
    s.set_visible(False)
ax.set_title("Cayley table of D4: entry = a∘b.  Blue = rotation, red = reflection", fontsize=10.5, pad=10, loc="center", y=-0.1)
save(fig, "ch56-cayley-table")
