"""Figures for Chapter 93 — Tensor Mathematics and Multilinear Algebra."""
import sys, pathlib; sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
from style import plt, np, PALETTE, INK, MUTED, new_figure, diagram, node, box, arrow, save
from matplotlib.patches import Polygon, Rectangle

BLUE, RED, GREEN, PURPLE, ORANGE, TEAL, GREY = PALETTE


# ---------------------------------------------------------------- figure 1: an image as a sum of outer products
def picture(n=32):
    y, x = np.mgrid[0:n, 0:n] / (n - 1)
    img = 0.15 + 0.5 * y                                          # sky gradient
    img = np.where((x > 0.2) & (x < 0.8) & (y > 0.5), 0.9, img)  # house body
    roof = (y > 0.2) & (y <= 0.5) & (np.abs(x - 0.5) < (y - 0.2) * 1.15)
    img = np.where(roof, 0.35, img)                               # triangular roof
    img = np.where((x > 0.43) & (x < 0.57) & (y > 0.72), 0.1, img)  # door
    img = np.where((x > 0.27) & (x < 0.37) & (y > 0.6) & (y < 0.7), 0.55, img)  # window
    return img


M = picture()
U, s, Vt = np.linalg.svd(M)
fig, axs = plt.subplots(2, 5, figsize=(9.6, 4.4), constrained_layout=True)
lim = max(abs(s[k] * np.outer(U[:, k], Vt[k])).max() for k in range(4))
for k in range(4):
    ax = axs[0, k]
    term = s[k] * np.outer(U[:, k], Vt[k])
    ax.imshow(term, cmap="RdBu_r", vmin=-lim, vmax=lim)
    ax.set_title(rf"$\sigma_{k+1}\,u_{k+1}v_{k+1}^\mathsf{{T}}$" + f"\nσ = {s[k]:.2f}", fontsize=9.5)
axs[0, 4].semilogy(np.arange(1, len(s) + 1), s, "o-", color=BLUE, ms=3, lw=1.2)
axs[0, 4].set_title("all singular values", fontsize=9.5)
axs[0, 4].set_xlabel("k", fontsize=9)
axs[0, 4].tick_params(labelsize=8)
for j, r in enumerate([1, 2, 4, 8]):
    Mr = (U[:, :r] * s[:r]) @ Vt[:r]
    err = np.linalg.norm(M - Mr) / np.linalg.norm(M)
    axs[1, j].imshow(Mr, cmap="gray", vmin=0, vmax=1)
    axs[1, j].set_title(f"sum of first {r}\nerror {100 * err:.1f}%", fontsize=9.5)
axs[1, 4].imshow(M, cmap="gray", vmin=0, vmax=1)
axs[1, 4].set_title(f"original\nrank {np.linalg.matrix_rank(M)}", fontsize=9.5)
for ax in list(axs[0, :4]) + list(axs[1]):
    ax.set_xticks([]); ax.set_yticks([]); ax.grid(False)
fig.suptitle("A 32×32 picture as a sum of rank-1 outer products (top: the terms, bottom: partial sums)",
             fontsize=11, fontweight="bold")
save(fig, "ch93-rank-one-image")
print("singular values", np.round(s[:10], 3), "rank", np.linalg.matrix_rank(M))
for r in (1, 2, 4, 8):
    Mr = (U[:, :r] * s[:r]) @ Vt[:r]
    print(r, np.linalg.norm(M - Mr) / np.linalg.norm(M))


# ---------------------------------------------------------------- figure 2: a 3x4x2 tensor and its unfoldings
I, J, K = 3, 4, 2
X = np.zeros((I, J, K), dtype=int)
for i in range(I):
    for j in range(J):
        for k in range(K):
            X[i, j, k] = 1 + i + I * j + I * J * k


def unfold(X, mode):
    return np.reshape(np.moveaxis(X, mode, 0), (X.shape[mode], -1), order="F")


fig, ax = plt.subplots(figsize=(9.6, 4.6), constrained_layout=True)
ax.set_xlim(0, 19.2); ax.set_ylim(0, 9.2); ax.set_aspect("equal"); ax.axis("off")
cell, dx, dy = 0.8, 1.1, 0.85                 # cell size and offset of the back slice
x0, y0 = 0.7, 2.9
slice_col = [BLUE, ORANGE]
for k in (1, 0):                               # draw the back slice first
    ox, oy = x0 + k * dx, y0 + k * dy
    for i in range(I):
        for j in range(J):
            cx, cy = ox + j * cell, oy + (I - 1 - i) * cell
            ax.add_patch(Rectangle((cx, cy), cell, cell, fc="white" if k == 0 else "#fff3e8",
                                   ec=slice_col[k], lw=1.3, zorder=3 - k))
            ax.text(cx + cell / 2, cy + cell / 2, str(X[i, j, k]), ha="center", va="center",
                    fontsize=9.5, color=INK, zorder=4 - k)
# depth edges
for (cx, cy) in [(x0, y0 + I * cell), (x0 + J * cell, y0 + I * cell), (x0 + J * cell, y0)]:
    ax.plot([cx, cx + dx], [cy, cy + dy], color=MUTED, lw=1)
ax.text(x0 + J * cell / 2, y0 - 0.35, "j →  (4 columns)", ha="center", fontsize=8.5, color=MUTED)
ax.text(x0 - 0.3, y0 + I * cell / 2, "i ↓", ha="right", va="center", fontsize=8.5, color=MUTED)
ax.text(x0 + J * cell + dx + 0.1, y0 + I * cell + dy / 2, "k ↗ (2 slices)", fontsize=8.5, color=MUTED)
ax.text(x0 + 2.0, 6.4, r"tensor $\mathcal{X}\in\mathbb{R}^{3\times4\times2}$", ha="center", fontsize=11)
ax.text(x0 + 1.6, 1.6, "front slice k = 1 (blue)\nback slice k = 2 (orange)", ha="center",
        fontsize=8.5, color=MUTED)


def table(ax, A, x, y, title, colfun):
    r, c = A.shape
    w = 0.5
    for a in range(r):
        for b in range(c):
            ax.add_patch(Rectangle((x + b * w, y + (r - 1 - a) * w), w, w, fc=colfun(a, b),
                                   ec="#9ca3af", lw=0.8))
            ax.text(x + b * w + w / 2, y + (r - 1 - a) * w + w / 2, str(A[a, b]),
                    ha="center", va="center", fontsize=7.5, color=INK)
    ax.text(x, y + r * w + 0.12, title, fontsize=9.5, color=INK)


light = {0: "white", 1: "#fff3e8"}
X1, X2, X3 = unfold(X, 0), unfold(X, 1), unfold(X, 2)
table(ax, X1, 6.4, 6.6, r"mode-1 unfolding $X_{(1)}$: 3×8 (rows = i)",
      lambda a, b: light[b // J])
table(ax, X2, 6.4, 1.7, r"mode-2 unfolding $X_{(2)}$: 4×6 (rows = j)",
      lambda a, b: light[b // I])
table(ax, X3, 12.4, 4.4, r"mode-3 unfolding $X_{(3)}$: 2×12 (rows = k)",
      lambda a, b: light[a])
for (tx, ty) in [(6.3, 7.3), (6.3, 2.7), (12.3, 4.9)]:
    ax.annotate("", xy=(tx, ty), xytext=(5.0, 4.6),
                arrowprops=dict(arrowstyle="-|>", color=MUTED, lw=1.2,
                                connectionstyle="arc3,rad=0.0"))
ax.text(14.0, 1.2, "Each unfolding lays the same 24 numbers out as a matrix.\n"
        "Their ranks (2, 2, 2 here) form the multilinear rank.", ha="center",
        fontsize=8.5, color=MUTED)
save(fig, "ch93-tensor-unfoldings")
print("unfolding ranks", [np.linalg.matrix_rank(u) for u in (X1, X2, X3)])
