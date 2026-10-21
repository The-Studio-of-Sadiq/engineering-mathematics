import sys, pathlib; sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
from style import plt, np, PALETTE, INK, MUTED, new_figure, diagram, node, box, arrow, save

# --- Figure 1: binary entropy and the capacity of a binary symmetric channel ----
p = np.linspace(1e-6, 1 - 1e-6, 600)
H = -p * np.log2(p) - (1 - p) * np.log2(1 - p)
fig, ax = new_figure(figsize=(7.2, 4.0))
ax.plot(p, H, color=PALETTE[0], label=r"entropy $H_b(p)$")
ax.plot(p, 1 - H, color=PALETTE[1], ls="--", label=r"BSC capacity $1-H_b(p)$")
for q in (0.11, 0.5):
    hq = -q * np.log2(q) - (1 - q) * np.log2(1 - q)
    ax.plot([q], [hq], "o", color=INK, zorder=5)
ax.annotate("fair coin: 1 bit,\nmaximal surprise", xy=(0.5, 1.0), xytext=(0.6, 0.4),
            arrowprops=dict(arrowstyle="->", color=INK), fontsize=9.5)
ax.annotate(r"$p=0.11$: half a bit", xy=(0.11, 0.5), xytext=(0.17, 0.26),
            arrowprops=dict(arrowstyle="->", color=INK), fontsize=9.5)
ax.annotate("useless channel", xy=(0.5, 0.0), xytext=(0.36, 0.12),
            arrowprops=dict(arrowstyle="->", color=PALETTE[1]), fontsize=9.5, color=PALETTE[1])
ax.set_xlabel(r"probability of heads (or of a bit flip) $p$")
ax.set_ylabel("bits")
ax.set_xlim(0, 1); ax.set_ylim(0, 1.08)
ax.set_title("One coin, two readings: uncertainty and capacity")
ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.16), ncol=2)
save(fig, "ch33-binary-entropy")

# --- Figure 2: a Huffman tree ----------------------------------------------------
fig, ax = diagram(figsize=(8.4, 4.8))
ax.set_ylim(-0.1, 6.0)
pos = {"root": (5.15, 5.4), "L": (2.65, 4.0), "R": (7.65, 4.0), "LL": (1.7, 2.6), "B": (3.6, 2.6),
       "RL": (6.7, 2.6), "A": (8.6, 2.6), "F": (1.0, 1.2), "E": (2.4, 1.2), "D": (6.0, 1.2), "C": (7.4, 1.2)}
label = {"root": "1.00", "L": "0.38", "R": "0.62", "LL": "0.18", "RL": "0.27",
         "A": "A\n0.35", "B": "B\n0.20", "C": "C\n0.15", "D": "D\n0.12", "E": "E\n0.10", "F": "F\n0.08"}
code = {"A": "11", "B": "01", "C": "101", "D": "100", "E": "001", "F": "000"}
edges = [("root", "L", 0), ("root", "R", 1), ("L", "LL", 0), ("L", "B", 1), ("R", "RL", 0), ("R", "A", 1),
         ("LL", "F", 0), ("LL", "E", 1), ("RL", "D", 0), ("RL", "C", 1)]
for a, b, bit in edges:
    (x0, y0), (x1, y1) = pos[a], pos[b]
    ax.plot([x0, x1], [y0, y1], color=MUTED, lw=1.6, zorder=1)
    ax.text((x0 + x1) / 2 + (-0.18 if bit == 0 else 0.18), (y0 + y1) / 2 + 0.08, str(bit),
            color=PALETTE[1], fontsize=11, weight="bold", ha="center")
for k, (x, y) in pos.items():
    leaf = k in code
    node(ax, x, y, label[k], color=PALETTE[2] if leaf else PALETTE[0], r=0.45, fontsize=9)
    if leaf:
        ax.text(x, y - 0.68, code[k], ha="center", fontsize=10.5, family="monospace", color=INK)
ax.text(5.0, -0.05, "average length 2.45 bits/symbol   vs   entropy 2.40 bits/symbol",
        ha="center", fontsize=10, color=INK)
save(fig, "ch33-huffman-tree")

# --- Figure 3: the entropy Venn diagram --------------------------------------------
fig, ax = diagram(figsize=(7.0, 3.8))
ax.set_ylim(0.6, 5.4)
from matplotlib.patches import Circle
ax.add_patch(Circle((4.0, 3.0), 2.0, fc=PALETTE[0], alpha=0.15, ec=PALETTE[0], lw=2))
ax.add_patch(Circle((6.0, 3.0), 2.0, fc=PALETTE[1], alpha=0.15, ec=PALETTE[1], lw=2))
ax.text(3.0, 3.0, r"$H(X\mid Y)$", ha="center", va="center", fontsize=12, color=INK)
ax.text(7.0, 3.0, r"$H(Y\mid X)$", ha="center", va="center", fontsize=12, color=INK)
ax.text(5.0, 3.0, r"$I(X;Y)$", ha="center", va="center", fontsize=12, color=INK, weight="bold")
ax.text(2.6, 5.15, r"$H(X)$", color=PALETTE[0], fontsize=12, ha="center")
ax.text(7.4, 5.15, r"$H(Y)$", color=PALETTE[1], fontsize=12, ha="center")
ax.text(5.0, 0.75, r"whole picture: $H(X,Y) = H(X)+H(Y)-I(X;Y)$", ha="center", fontsize=10.5, color=MUTED)
save(fig, "ch33-entropy-venn")
