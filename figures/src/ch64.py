import sys, pathlib; sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
from style import plt, np, PALETTE, INK, MUTED, new_figure, diagram, node, box, arrow, save
from matplotlib.patches import Rectangle, FancyBboxPatch
from math import gcd

# ---------------------------------------------------------------- 1. Cantor's diagonal argument
rng = np.random.default_rng(64)
R, C = 8, 10
bits = rng.integers(0, 2, size=(R, C))
fig, ax = plt.subplots(figsize=(8.6, 5.6), constrained_layout=True)
ax.set_xlim(-2.6, C + 0.4); ax.set_ylim(-2.3, R + 1.0); ax.set_aspect("equal"); ax.axis("off")
for i in range(R):
    ax.text(-0.25, R - 1 - i + 0.5, f"$s_{{{i+1}}}$", ha="right", va="center", fontsize=12, color=INK)
    for k in range(C):
        diag = (i == k)
        fc = "#fde68a" if diag else "white"
        ax.add_patch(Rectangle((k, R - 1 - i), 1, 1, fc=fc, ec="#9ca3af", lw=1.0))
        ax.text(k + 0.5, R - 1 - i + 0.5, str(bits[i, k]), ha="center", va="center", fontsize=12,
                color=PALETTE[1] if diag else INK, weight="bold" if diag else "normal")
    ax.text(C + 0.2, R - 1 - i + 0.5, r"$\cdots$", ha="left", va="center", fontsize=12, color=MUTED)
for k in range(C):
    ax.text(k + 0.5, R + 0.3, str(k + 1), ha="center", va="bottom", fontsize=9, color=MUTED)
ax.text(C / 2, R + 0.75, "digit position", ha="center", fontsize=9.5, color=MUTED)
ax.text(-1.4, R / 2, "the list\n(any list)", ha="center", va="center", fontsize=9.5, color=MUTED, rotation=90)
ax.plot([0, R], [R, 0], color=PALETTE[1], lw=1.2, ls="--", alpha=0.6)
# the new sequence
d = [1 - bits[k, k] if k < R else int(rng.integers(0, 2)) for k in range(C)]
for k in range(C):
    fc = "#dcfce7" if k < R else "#f3f4f6"
    ax.add_patch(Rectangle((k, -1.8), 1, 1, fc=fc, ec=PALETTE[2], lw=1.6))
    ax.text(k + 0.5, -1.3, str(d[k]), ha="center", va="center", fontsize=12, color=PALETTE[2], weight="bold")
ax.text(-0.25, -1.3, "$d$", ha="right", va="center", fontsize=13, color=PALETTE[2])
ax.text(C / 2, -0.45, r"flip every diagonal digit:  $d_k = 1 - s_k(k)$,  so $d \neq s_k$ for every $k$",
        ha="center", va="center", fontsize=10, color=PALETTE[2])
save(fig, "ch64-cantor-diagonal")

# ---------------------------------------------------------------- 2. Hilbert's hotel
fig, ax = diagram(figsize=(10.0, 4.6))
ax.set_xlim(0, 13); ax.set_ylim(0, 6)
rooms = 9
for row, (y, title) in enumerate([(4.3, "One new guest:  n → n + 1"),
                                   (1.4, "Infinitely many new guests:  n → 2n")]):
    ax.text(0.2, y + 1.15, title, fontsize=10.5, color=INK, weight="bold")
    for k in range(rooms):
        x = 0.6 + 1.3 * k
        ax.add_patch(FancyBboxPatch((x, y - 0.35), 0.9, 0.7, boxstyle="round,pad=0.02,rounding_size=0.08",
                                    fc="#eef2ff", ec=PALETTE[0], lw=1.4))
        ax.text(x + 0.45, y, str(k + 1), ha="center", va="center", fontsize=10, color=INK)
        tgt = (k + 1) if row == 0 else (2 * k + 1)
        if tgt < rooms:
            ax.annotate("", xy=(0.6 + 1.3 * tgt + 0.45, y + 0.38), xytext=(x + 0.45, y + 0.38),
                        arrowprops=dict(arrowstyle="-|>", color=PALETTE[0], lw=1.1,
                                        connectionstyle="arc3,rad=-0.5"))
    ax.text(0.6 + 1.3 * rooms + 0.1, y, r"$\cdots$", fontsize=14, color=MUTED, va="center")
    if row == 0:
        ax.text(1.05, y - 0.75, "free!", ha="center", fontsize=9.5, color=PALETTE[2])
    else:
        for k in range(0, rooms, 2):
            ax.text(0.6 + 1.3 * k + 0.45, y - 0.75, "free", ha="center", fontsize=9, color=PALETTE[2])
save(fig, "ch64-hilbert-hotel")

# ---------------------------------------------------------------- 3. zig-zag enumeration of the positive rationals
N = 6
fig, ax = plt.subplots(figsize=(7.2, 6.2), constrained_layout=True)
ax.set_xlim(0.3, N + 0.8); ax.set_ylim(N + 0.7, 0.3); ax.set_aspect("equal"); ax.axis("off")
order = []
for s in range(2, 2 * N + 1):
    diag = [(p, s - p) for p in range(1, s) if p <= N and s - p <= N]
    if s % 2 == 1:
        diag = diag[::-1]
    order.extend(diag)
order = [pq for pq in order if pq[0] + pq[1] <= N + 1]
count = 0
labels = {}
for p, q in order:
    if gcd(p, q) == 1:
        count += 1
        labels[(p, q)] = count
for q in range(1, N + 1):
    for p in range(1, N + 1):
        new = gcd(p, q) == 1
        col = INK if new else "#c0c4cc"
        ax.text(p, q, f"{p}/{q}", ha="center", va="center", fontsize=12, color=col,
                bbox=dict(boxstyle="round,pad=0.25", fc="#eef2ff" if (p, q) in labels else "white",
                          ec=PALETTE[0] if (p, q) in labels else "#e5e7eb", lw=1.0))
        if (p, q) in labels:
            ax.text(p - 0.3, q - 0.27, str(labels[(p, q)]), ha="right", fontsize=8.5, color=PALETTE[1], weight="bold")
for (a, b), (c, d) in zip(order[:-1], order[1:]):
    ax.annotate("", xy=(c, d), xytext=(a, b),
                arrowprops=dict(arrowstyle="-|>", color=PALETTE[4], lw=1.2, shrinkA=13, shrinkB=13, alpha=0.85))
for p in range(1, N + 1):
    ax.text(p, 0.45, f"p = {p}", ha="center", fontsize=9, color=MUTED)
for q in range(1, N + 1):
    ax.text(0.35, q, f"q = {q}", ha="left", va="center", fontsize=9, color=MUTED, rotation=0)
ax.text((N + 1) / 2 + 0.3, N + 0.55, "grey entries (2/2, 2/4, ...) repeat an earlier value and are skipped",
        ha="center", fontsize=9, color=MUTED)
save(fig, "ch64-rationals-zigzag")
