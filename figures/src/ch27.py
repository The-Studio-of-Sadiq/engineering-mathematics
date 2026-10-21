import sys, pathlib; sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
from style import plt, np, PALETTE, INK, MUTED, new_figure, diagram, node, box, arrow, save
import heapq

# ---------------------------------------------------------------- weighted road network + Dijkstra
pos = {"A": (0.8, 3.0), "B": (3.2, 5.0), "C": (3.0, 1.2), "D": (5.6, 3.4),
       "E": (6.2, 0.9), "F": (8.2, 5.1), "G": (9.2, 2.2)}
edges = [("A", "B", 4), ("A", "C", 2), ("C", "B", 1), ("B", "D", 5), ("C", "D", 8),
         ("C", "E", 10), ("D", "E", 2), ("D", "F", 6), ("E", "G", 3), ("F", "G", 4)]
adj = {v: [] for v in pos}
for u, v, w in edges:
    adj[u].append((v, w)); adj[v].append((u, w))
dist, prev, pq = {"A": 0}, {}, [(0, "A")]
while pq:
    d, u = heapq.heappop(pq)
    if d > dist[u]:
        continue
    for v, w in adj[u]:
        if d + w < dist.get(v, 1e9):
            dist[v], prev[v] = d + w, u
            heapq.heappush(pq, (d + w, v))
path, v = ["G"], "G"
while v != "A":
    v = prev[v]; path.append(v)
on_path = {frozenset(e) for e in zip(path, path[1:])}

fig, ax = diagram(figsize=(8.0, 4.8))
ax.set_xlim(0, 10); ax.set_ylim(-0.3, 6.0)
for u, v, w in edges:
    (x0, y0), (x1, y1) = pos[u], pos[v]
    hot = frozenset((u, v)) in on_path
    ax.plot([x0, x1], [y0, y1], color=PALETTE[1] if hot else "#c7cdd6", lw=4 if hot else 2, zorder=1,
            solid_capstyle="round")
    mx, my = (x0 + x1) / 2, (y0 + y1) / 2
    ax.text(mx, my, str(w), ha="center", va="center", fontsize=10, color=PALETTE[1] if hot else INK,
            bbox=dict(boxstyle="round,pad=0.18", fc="white", ec="none"), zorder=2)
for k, (x, y) in pos.items():
    node(ax, x, y, k, color=PALETTE[1] if k in path else PALETTE[0], fc="#fff1f1" if k in path else "white")
    lx, ly = {"B": (0, 0.62), "F": (0, 0.62), "E": (0.8, -0.3), "D": (0.75, -0.45)}.get(k, (0, -0.68))
    ax.text(x + lx, y + ly, f"d={dist[k]}", ha="center", fontsize=8.5, color=MUTED)
ax.text(5.0, -0.2, "shortest A to G: " + " - ".join(reversed(path)) + f"   (length {dist['G']})",
        ha="center", fontsize=10, color=PALETTE[1])
save(fig, "ch27-shortest-path")
print("path", list(reversed(path)), dist)

# ---------------------------------------------------------------- bipartite matching: greedy vs augmented
D = ["D1", "D2", "D3", "D4", "D5"]
R = ["R1", "R2", "R3", "R4", "R5"]
E = [("D1", "R1"), ("D1", "R2"), ("D2", "R1"), ("D3", "R2"), ("D3", "R3"),
     ("D4", "R3"), ("D4", "R5"), ("D5", "R4"), ("D5", "R5")]
greedy = [("D1", "R1"), ("D3", "R2"), ("D4", "R3"), ("D5", "R4")]
augpath = [("D2", "R1"), ("D1", "R1"), ("D1", "R2"), ("D3", "R2"), ("D3", "R3"), ("D4", "R3"), ("D4", "R5")]
best = [("D1", "R2"), ("D2", "R1"), ("D3", "R3"), ("D4", "R5"), ("D5", "R4")]

fig, axes = plt.subplots(1, 2, figsize=(9.6, 4.6), constrained_layout=True)
for ax, title, M, P in [(axes[0], "Greedy: 4 pairs, D2 stranded", greedy, augpath),
                        (axes[1], "After flipping the path: 5 pairs", best, None)]:
    ax.set_xlim(0, 6); ax.set_ylim(-0.3, 5.9); ax.set_aspect("equal"); ax.axis("off")
    ax.set_title(title, fontsize=11, color=INK)
    py = {n: 5 - i for i, n in enumerate(D)}
    py.update({n: 5 - i for i, n in enumerate(R)})
    px = {n: 1.2 for n in D}; px.update({n: 4.8 for n in R})
    for u, v in E:
        inM = (u, v) in M
        inP = P is not None and (u, v) in P
        col = PALETTE[2] if inM else ("#c7cdd6")
        ax.plot([px[u], px[v]], [py[u], py[v]], color=col, lw=4 if inM else 1.6,
                ls="-" if inM or not inP else "--", zorder=1)
        if inP and not inM:
            ax.plot([px[u], px[v]], [py[u], py[v]], color=PALETTE[4], lw=2.4, ls="--", zorder=1)
        if inP and inM:
            ax.plot([px[u], px[v]], [py[u], py[v]], color=PALETTE[4], lw=1.2, zorder=2)
    for n in D:
        stranded = (M is greedy and n == "D2")
        node(ax, px[n], py[n], n, color=PALETTE[1] if stranded else PALETTE[0], r=0.38, fontsize=9)
    for n in R:
        free = (M is greedy and n == "R5")
        node(ax, px[n], py[n], n, color=PALETTE[1] if free else PALETTE[3], r=0.38, fontsize=9)
    ax.text(1.2, -0.25, "drivers", ha="center", color=PALETTE[0], fontsize=9.5)
    ax.text(4.8, -0.25, "riders", ha="center", color=PALETTE[3], fontsize=9.5)
axes[0].text(3.0, 5.62, "orange: augmenting path D2-R1-D1-R2-D3-R3-D4-R5", ha="center", fontsize=8.3,
             color=PALETTE[4])
save(fig, "ch27-matching")

# ---------------------------------------------------------------- spectral communities
rng = np.random.default_rng(11)
n1 = n2 = 12
N = n1 + n2
grp = np.array([0] * n1 + [1] * n2)
P = np.where(grp[:, None] == grp[None, :], 0.45, 0.05)
A = (rng.random((N, N)) < P).astype(float)
A = np.triu(A, 1); A = A + A.T
L = np.diag(A.sum(1)) - A
vals, vecs = np.linalg.eigh(L)
fied = vecs[:, 1]
side = fied > 0
# layout: x from Fiedler value, y spread
k = np.arange(N) % n1
ang = 2 * np.pi * k / n1 + 0.3
rad = np.where(k % 2 == 0, 2.0, 1.15)
xy = np.c_[np.where(grp == 0, 2.5, 7.5) + rad * np.cos(ang), 3.0 + rad * np.sin(ang)]
fig, (ax, ax2) = plt.subplots(1, 2, figsize=(10.0, 4.2), constrained_layout=True,
                              gridspec_kw=dict(width_ratios=[1.35, 1]))
ax.set_xlim(0, 10); ax.set_ylim(0, 6); ax.set_aspect("equal"); ax.axis("off")
cross = 0
for i in range(N):
    for j in range(i + 1, N):
        if A[i, j]:
            between = side[i] != side[j]
            cross += between
            ax.plot(*zip(xy[i], xy[j]), color=PALETTE[1] if between else "#c7cdd6",
                    lw=1.6 if between else 1.0, zorder=1)
for i in range(N):
    c = PALETTE[0] if side[i] else PALETTE[4]
    ax.add_patch(plt.Circle(xy[i], 0.2, fc=c, ec="white", lw=1, zorder=3))
ax.set_title(f"Sign of the Fiedler vector splits the graph ({cross} red edges cut)", fontsize=10.5)
order = np.argsort(fied)
ax2.bar(range(N), fied[order], color=[PALETTE[0] if s else PALETTE[4] for s in side[order]])
ax2.axhline(0, color=INK, lw=0.8)
ax2.set_xlabel("vertices, sorted by entry"); ax2.set_xticks([])
ax2.set_ylabel("Fiedler vector entry")
ax2.set_title(fr"$\lambda_2 = {vals[1]:.2f}$ (second-smallest of $\mathbf{{L}}$)", fontsize=10.5)
save(fig, "ch27-spectral-communities")
print("lambda", vals[:3], "cut", cross, "correct", max((side == grp).mean(), (side != grp).mean()))
