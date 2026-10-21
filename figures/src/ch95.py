"""Figures for Chapter 95 — Representation Theory: Making Abstract Symmetry Concrete."""
import sys, pathlib; sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
from style import plt, np, PALETTE, INK, MUTED, new_figure, diagram, node, box, arrow, save

BLUE, RED, GREEN, PURPLE, ORANGE, TEAL, GREY = PALETTE

# ---------------------------------------------------------------- figure 1: character table of S3
fig, ax = plt.subplots(figsize=(7.4, 3.6), constrained_layout=True)
ax.set_xlim(0, 10); ax.set_ylim(0, 5); ax.axis("off")
cols = ["$e$", "$(1\\,2)$", "$(1\\,2\\,3)$"]
sizes = ["1 element", "3 elements", "2 elements"]
rows = [("trivial  $\\chi_1$", [1, 1, 1], BLUE),
        ("sign  $\\chi_2$", [1, -1, 1], GREEN),
        ("standard  $\\chi_3$", [2, 0, -1], RED)]
x0, cw, rh = 3.2, 1.9, 0.72
top = 4.3
ax.text(1.6, top, "irreducible", ha="center", va="center", fontsize=10, color=MUTED)
for k, (c, s) in enumerate(zip(cols, sizes)):
    cx = x0 + cw * (k + 0.5)
    ax.text(cx, top + 0.12, c, ha="center", va="center", fontsize=13, color=INK)
    ax.text(cx, top - 0.3, s, ha="center", va="center", fontsize=8.5, color=MUTED)
ax.plot([0.2, x0 + 3 * cw], [top - 0.6, top - 0.6], color=INK, lw=1.4)
ax.plot([x0, x0], [top - 0.6 - 3 * rh, top + 0.5], color=INK, lw=1.0)
for r, (name, vals, col) in enumerate(rows):
    cy = top - 0.6 - rh * (r + 0.5)
    ax.text(1.6, cy, name, ha="center", va="center", fontsize=11, color=col)
    for k, v in enumerate(vals):
        cx = x0 + cw * (k + 0.5)
        ax.add_patch(plt.Rectangle((cx - 0.55, cy - 0.27), 1.1, 0.54, fc=col, alpha=0.10 + 0.10 * abs(v), ec="none"))
        ax.text(cx, cy, f"{v:d}", ha="center", va="center", fontsize=14, color=INK)
ax.text(x0 + 1.5 * cw, 0.55,
        r"rows are orthogonal: $\langle\chi_2,\chi_3\rangle=\frac{1}{6}\,(1\cdot2\cdot1+3\cdot(-1)\cdot0+2\cdot1\cdot(-1))=0$"
        "\n" r"first column: $1^2+1^2+2^2=6=|S_3|$", ha="center", va="center", fontsize=9.5, color=INK)
ax.text(9.95, 2.6, "class\nsizes\nweight\nthe sum", ha="right", va="center", fontsize=8, color=MUTED)
save(fig, "ch95-s3-character-table")

# ---------------------------------------------------------------- figure 2: D4 on a square, and its pieces
fig, ax = plt.subplots(figsize=(9.6, 3.6), constrained_layout=True)
ax.set_xlim(0, 19.2); ax.set_ylim(0, 7.2); ax.set_aspect("equal"); ax.axis("off")
# the square with its symmetry axes
c = np.array([2.9, 3.6]); h = 1.6
corners = [c + h * np.array(v) for v in [(1, 1), (-1, 1), (-1, -1), (1, -1)]]
sq = np.array(corners + [corners[0]])
ax.plot(sq[:, 0], sq[:, 1], color=INK, lw=2)
ax.plot([c[0] - 2.4, c[0] + 2.4], [c[1] - 2.4, c[1] + 2.4], color=PURPLE, ls="--", lw=1.2)
ax.plot([c[0] - 2.4, c[0] + 2.4], [c[1] + 2.4, c[1] - 2.4], color=PURPLE, ls="--", lw=1.2)
ax.plot([c[0] - 2.4, c[0] + 2.4], [c[1], c[1]], color=ORANGE, ls=":", lw=1.6)
ax.plot([c[0], c[0]], [c[1] - 2.4, c[1] + 2.4], color=ORANGE, ls=":", lw=1.6)
for i, p in enumerate(corners):
    node(ax, p[0], p[1], str(i + 1), color=BLUE, r=0.36, fontsize=11)
ax.annotate("", xy=(c[0] + 0.0, c[1] + 0.75), xytext=(c[0] + 0.75, c[1] + 0.0),
            arrowprops=dict(arrowstyle="-|>", color=GREEN, lw=1.6, connectionstyle="arc3,rad=0.45"))
ax.text(c[0] + 0.62, c[1] + 0.62, "$r$", color=GREEN, fontsize=12)
ax.text(c[0], 0.45, "8 symmetries: 4 rotations,\n2 diagonal (dashed) + 2 edge (dotted) flips", ha="center",
        fontsize=8.5, color=MUTED)
ax.text(6.35, 3.6, r"$\mathbb{C}^4=$", fontsize=15, ha="center", va="center", color=INK)


def mini(cx, cy, vals, title, col):
    s = 0.85
    pts = [(cx + s, cy + s), (cx - s, cy + s), (cx - s, cy - s), (cx + s, cy - s)]
    P = np.array(pts + [pts[0]])
    ax.plot(P[:, 0], P[:, 1], color=MUTED, lw=1.0)
    for (px, py), v in zip(pts, vals):
        fc = {1: "#dbeafe", -1: "#fee2e2", 0: "white"}[v]
        ax.add_patch(plt.Circle((px, py), 0.33, fc=fc, ec=col, lw=1.4, zorder=3))
        ax.text(px, py, {1: "+1", -1: "−1", 0: "0"}[v], ha="center", va="center", fontsize=8.5, zorder=4, color=INK)
    ax.text(cx, cy - 1.75, title, ha="center", va="center", fontsize=9, color=col)


mini(8.4, 3.6, [1, 1, 1, 1], "trivial (dim 1)\nall moves: $+1$", BLUE)
ax.text(9.95, 3.6, r"$\oplus$", fontsize=15, ha="center", va="center")
mini(11.5, 3.6, [1, -1, 1, -1], "sign-like (dim 1)\n$r\\mapsto-1$, diag flip $\\mapsto+1$", GREEN)
ax.text(13.05, 3.6, r"$\oplus$", fontsize=15, ha="center", va="center")
mini(14.6, 3.6, [1, 0, -1, 0], "", RED)
mini(17.4, 3.6, [0, 1, 0, -1], "", RED)
ax.text(16.0, 3.6, ",", fontsize=15, ha="center", va="center", color=RED)
ax.text(16.0, 1.6, "2-dim irreducible:\nrotation by $r$ mixes the two", ha="center", fontsize=9, color=RED)
ax.text(9.6, 6.75, r"character of the corner permutations: $\chi(e)=4,\ \chi(r^2)=\chi(r)=0,\ \chi(\mathrm{edge})=0,\ \chi(\mathrm{diag})=2$",
        ha="center", fontsize=9.5, color=INK)
save(fig, "ch95-d4-square")

# ---------------------------------------------------------------- figure 3: Fourier on Z_n diagonalises convolution
n = 8
h_ = np.array([3, 1, 0, 0, 0, 0, 0, 1.0]) / 5  # symmetric blur kernel
C = np.array([[h_[(i - k) % n] for k in range(n)] for i in range(n)])
F = np.exp(-2j * np.pi * np.outer(np.arange(n), np.arange(n)) / n) / np.sqrt(n)
D = F @ C @ F.conj().T
fig, axs = new_figure(1, 3, figsize=(9.6, 3.3))
im = axs[0].imshow(C, cmap="Blues", vmin=0, vmax=0.6)
axs[0].set_title("circular convolution $\\mathbf{C}$\n(commutes with every shift)", fontsize=10)
axs[1].imshow(np.abs(D), cmap="Reds", vmin=0, vmax=1)
axs[1].set_title("$\\mathbf{F}\\,\\mathbf{C}\\,\\mathbf{F}^{\\mathsf{H}}$: diagonal\n(entries $=\\hat h[k]$)", fontsize=10)
for a in axs[:2]:
    a.set_xticks(range(n)); a.set_yticks(range(n)); a.grid(False)
    a.tick_params(labelsize=7)
m = np.arange(n); mm = np.linspace(0, n, 400)
for k, col in zip([0, 1, 2], [BLUE, GREEN, RED]):
    axs[2].plot(mm, np.cos(2 * np.pi * k * mm / n) + 2.6 * (2 - k), color=col, lw=1, alpha=0.35)
    axs[2].plot(m, np.cos(2 * np.pi * k * m / n) + 2.6 * (2 - k), "o", color=col, ms=5)
    axs[2].text(8.25, 2.6 * (2 - k) + 1, f"$k={k}$", va="center", color=col, fontsize=9)
axs[2].set_yticks([]); axs[2].set_xticks(range(n)); axs[2].set_xlim(-0.4, 9.3)
axs[2].set_xlabel("$m\\in\\mathbb{Z}_8$")
axs[2].set_title("characters $\\operatorname{Re}\\,e^{2\\pi jkm/8}$:\nthe irreducible pieces", fontsize=10)
save(fig, "ch95-zn-fourier")
