import sys, pathlib; sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
from style import plt, np, PALETTE, INK, MUTED, new_figure, diagram, node, box, arrow, save
from matplotlib.colors import ListedColormap


def eca(rule, width=161, steps=80):
    """Elementary cellular automaton from a single live cell, periodic boundary."""
    table = np.array([(rule >> k) & 1 for k in range(8)], dtype=np.uint8)
    x = np.zeros(width, dtype=np.uint8)
    x[width // 2] = 1
    rows = [x]
    for _ in range(steps - 1):
        idx = 4 * np.roll(x, 1) + 2 * x + np.roll(x, -1)
        x = table[idx]
        rows.append(x)
    return np.array(rows)


# --- Figure 1: space-time diagrams of rules 30, 90, 110 --------------------------------
fig, axes = new_figure(1, 3, figsize=(10.5, 4.0))
specs = [(30, "Rule 30: chaos", PALETTE[1]),
         (90, "Rule 90: Sierpinski fractal", PALETTE[0]),
         (110, "Rule 110: gliders that compute", PALETTE[2])]
for ax, (r, title, col) in zip(axes, specs):
    g = eca(r)
    if r == 110:   # rule 110 grows to the left only: shift the seed to the right edge
        g = eca(r, width=161, steps=80)
        g = np.roll(g, 60, axis=1)
    cmap = ListedColormap(["white", col])
    ax.imshow(g, cmap=cmap, interpolation="nearest", aspect="auto")
    ax.set_title(title, fontsize=11)
    ax.set_xlabel("cell position $i$")
    ax.grid(False)
    ax.set_xticks([])
    ax.set_yticks([0, 20, 40, 60, 79])
axes[0].set_ylabel("time step $t$  (downwards)")
save(fig, "ch86-elementary-rules")


# --- Figure 2: Game of Life glider ---------------------------------------------------
def life_step(G):
    n = sum(np.roll(np.roll(G, a, 0), b, 1) for a in (-1, 0, 1) for b in (-1, 0, 1)) - G
    return ((n == 3) | ((G == 1) & (n == 2))).astype(np.uint8)


G = np.zeros((8, 8), dtype=np.uint8)
for (r, c) in [(1, 2), (2, 3), (3, 1), (3, 2), (3, 3)]:
    G[r, c] = 1
fig, axes = new_figure(1, 5, figsize=(10.5, 2.6))
for t, ax in enumerate(axes):
    ax.imshow(G, cmap=ListedColormap(["white", PALETTE[3]]), interpolation="nearest")
    ax.set_xticks(np.arange(-0.5, 8, 1), minor=True)
    ax.set_yticks(np.arange(-0.5, 8, 1), minor=True)
    ax.grid(which="minor", color="#d1d5db", lw=0.8)
    ax.grid(which="major", visible=False)
    ax.tick_params(which="both", length=0, labelbottom=False, labelleft=False)
    for s in ax.spines.values():
        s.set_visible(True)
        s.set_color(MUTED)
    ax.set_title(f"generation {t}", fontsize=10.5)
    G = life_step(G)
fig.suptitle("After 4 generations the same shape reappears, shifted one cell down and one cell right",
             fontsize=10.5, color=INK)
save(fig, "ch86-glider")
