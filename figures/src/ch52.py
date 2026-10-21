"""Figures for Chapter 52 — Dynamic Programming and Bellman's Principle."""
import sys, pathlib; sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
from style import plt, np, PALETTE, INK, MUTED, new_figure, diagram, node, box, arrow, save

BLUE, RED, GREEN, PURPLE, ORANGE, TEAL, GREY = PALETTE

# ---------------------------------------------------------------- figure 1: staged shortest path, solved backwards
layers = [["S"], ["A1", "A2", "A3"], ["B1", "B2", "B3"], ["C1", "C2"], ["T"]]
w = {("S", "A1"): 2, ("S", "A2"): 4, ("S", "A3"): 4,
     ("A1", "B1"): 7, ("A1", "B2"): 6, ("A1", "B3"): 8,
     ("A2", "B1"): 3, ("A2", "B2"): 2, ("A2", "B3"): 4,
     ("A3", "B1"): 4, ("A3", "B2"): 1, ("A3", "B3"): 5,
     ("B1", "C1"): 1, ("B1", "C2"): 4, ("B2", "C1"): 6, ("B2", "C2"): 3,
     ("B3", "C1"): 3, ("B3", "C2"): 3, ("C1", "T"): 3, ("C2", "T"): 4}
V = {"T": 0}
best = {}
for layer in reversed(layers[:-1]):
    for u in layer:
        opts = [(w[(u, v)] + V[v], v) for v in [v for (a, v) in w if a == u]]
        V[u], best[u] = min(opts)
path = ["S"]
while path[-1] != "T":
    path.append(best[path[-1]])
greedy = ["S"]
while greedy[-1] != "T":
    u = greedy[-1]
    greedy.append(min([v for (a, v) in w if a == u], key=lambda v: w[(u, v)]))
gcost = sum(w[(greedy[i], greedy[i + 1])] for i in range(len(greedy) - 1))
if __name__ == "__main__":
    print("cost-to-go", V, "optimal path", path, "greedy", greedy, gcost)

pos = {}
for i, layer in enumerate(layers):
    ys = np.linspace(4.6, 1.4, len(layer)) if len(layer) > 1 else [3.0]
    for name, y in zip(layer, ys):
        pos[name] = (0.9 + 2.05 * i, y)

fig, ax = diagram(figsize=(9.4, 5.0))
ax.set_xlim(0, 10.2); ax.set_ylim(0.3, 6.0)
on_path = set(zip(path[:-1], path[1:]))
on_greedy = set(zip(greedy[:-1], greedy[1:]))
for (u, v), c in w.items():
    (x0, y0), (x1, y1) = pos[u], pos[v]
    if (u, v) in on_path:
        col, lw = GREEN, 3.0
    elif (u, v) in on_greedy:
        col, lw = ORANGE, 2.2
    else:
        col, lw = "#cbd5e1", 1.1
    ax.plot([x0, x1], [y0, y1], color=col, lw=lw, zorder=1, ls="--" if (u, v) in on_greedy and (u, v) not in on_path else "-")
    t = 0.26
    mx, my = x0 + t * (x1 - x0), y0 + t * (y1 - y0)
    ax.text(mx, my, str(c), fontsize=8.5, ha="center", va="center", color=INK, zorder=4,
            bbox=dict(boxstyle="round,pad=0.12", fc="white", ec="none"))
for name, (x, y) in pos.items():
    col = GREEN if name in path else BLUE
    node(ax, x, y, name, color=col, r=0.34, fontsize=9.5)
    ax.text(x, y + 0.5, f"V={V[name]}", ha="center", va="bottom", fontsize=9, color=PURPLE, weight="bold")
for i, lab in enumerate(["start", "stage 1", "stage 2", "stage 3", "goal"]):
    ax.text(0.9 + 2.05 * i, 0.55, lab, ha="center", fontsize=9, color=MUTED)
ax.plot([], [], color=GREEN, lw=3, label=f"optimal path, cost {V['S']}")
ax.plot([], [], color=ORANGE, lw=2.2, ls="--", label=f"greedy path, cost {gcost}")
ax.legend(loc="upper center", ncol=2, fontsize=9.5, bbox_to_anchor=(0.5, 1.02))
save(fig, "ch52-staged-shortest-path")

# ---------------------------------------------------------------- figure 2: grid-world values and policy
ROWS, COLS = 5, 6
walls = {(1, 1), (1, 2), (3, 3)}
goal, pit = (0, 5), (1, 5)
step, gamma, slip = -0.04, 0.95, 0.2
moves = {"U": (-1, 0), "D": (1, 0), "L": (0, -1), "R": (0, 1)}
perp = {"U": "LR", "D": "LR", "L": "UD", "R": "UD"}
cells = [(r, c) for r in range(ROWS) for c in range(COLS) if (r, c) not in walls]


def nxt(s, a):
    r, c = s[0] + moves[a][0], s[1] + moves[a][1]
    return (r, c) if 0 <= r < ROWS and 0 <= c < COLS and (r, c) not in walls else s


def q(s, a, Vg):
    out = (1 - slip) * Vg[nxt(s, a)]
    for b in perp[a]:
        out += slip / 2 * Vg[nxt(s, b)]
    return step + gamma * out


Vg = {s: 0.0 for s in cells}
for it in range(500):
    new = {}
    for s in cells:
        if s == goal:
            new[s] = 1.0
        elif s == pit:
            new[s] = -1.0
        else:
            new[s] = max(q(s, a, Vg) for a in moves)
    delta = max(abs(new[s] - Vg[s]) for s in cells)
    Vg = new
    if delta < 1e-10:
        break
pol = {s: max(moves, key=lambda a: q(s, a, Vg)) for s in cells if s not in (goal, pit)}
if __name__ == "__main__":
    print("iterations", it, {s: round(v, 3) for s, v in Vg.items()})

fig, ax = new_figure(figsize=(7.4, 5.4))
M = np.full((ROWS, COLS), np.nan)
for s, v in Vg.items():
    M[s] = v
cmap = plt.get_cmap("RdYlGn").copy(); cmap.set_bad("#374151")
im = ax.imshow(M, cmap=cmap, vmin=-1, vmax=1)
for (r, c) in walls:
    ax.text(c, r, "wall", ha="center", va="center", color="white", fontsize=9)
for s, v in Vg.items():
    r, c = s
    if s == goal:
        ax.text(c, r, "GOAL\n+1", ha="center", va="center", fontsize=10, weight="bold", color=INK)
        continue
    if s == pit:
        ax.text(c, r, "PIT\n−1", ha="center", va="center", fontsize=10, weight="bold", color="white")
        continue
    dr, dc = moves[pol[s]]
    ax.annotate("", xy=(c + 0.28 * dc, r + 0.28 * dr - 0.08), xytext=(c - 0.28 * dc, r - 0.28 * dr - 0.08),
                arrowprops=dict(arrowstyle="-|>", color=INK, lw=1.8, mutation_scale=16))
    ax.text(c, r + 0.33, f"{v:.2f}", ha="center", va="center", fontsize=8.5, color=INK)
ax.set_xticks(np.arange(-0.5, COLS, 1), minor=True); ax.set_yticks(np.arange(-0.5, ROWS, 1), minor=True)
ax.grid(which="minor", color="white", lw=2); ax.grid(which="major", visible=False)
ax.set_xticks([]); ax.set_yticks([])
ax.tick_params(which="minor", length=0)
for sp in ax.spines.values():
    sp.set_visible(False)
cb = fig.colorbar(im, ax=ax, shrink=0.8, pad=0.02)
cb.set_label(r"optimal value $V^*(s)$")
ax.set_title(r"value iteration: $\gamma=0.95$, step reward $-0.04$, slip 20%", fontsize=10.5)
save(fig, "ch52-gridworld-values")
