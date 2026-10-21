"""Figures for Chapter 77 — Cryptography: Mathematics Against Adversaries."""
import sys, pathlib; sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
from style import plt, np, PALETTE, INK, MUTED, new_figure, diagram, node, box, arrow, save

BLUE, RED, GREEN, PURPLE, ORANGE, TEAL, GREY = PALETTE

# ---------------------------------------------------------------- figure 1: symmetric vs public key
fig, ax = diagram(figsize=(9.6, 6.9))
ax.set_xlim(0, 12); ax.set_ylim(0, 8.7)

ax.text(0.2, 8.45, "(a) symmetric: one shared secret key k", fontsize=11, weight="bold", color=INK)
y = 6.5
node(ax, 1.0, y, "Alice", color=BLUE, r=0.55, fontsize=9.5)
box(ax, 3.4, y, "Enc$_k(m)$", w=1.7, h=0.8, color=BLUE)
box(ax, 8.6, y, "Dec$_k(c)$", w=1.7, h=0.8, color=GREEN)
node(ax, 11.0, y, "Bob", color=GREEN, r=0.55, fontsize=9.5)
arrow(ax, 1.0, y, 3.4, y, text="m", shrink=30)
arrow(ax, 3.4, y, 8.6, y, text="ciphertext c (public channel)", shrink=45, color=INK)
arrow(ax, 8.6, y, 11.0, y, text="m", shrink=40)
node(ax, 6.0, 5.0, "Eve", color=RED, r=0.5, fontsize=9.5)
ax.plot([6.0, 6.0], [5.5, 6.4], color=RED, lw=1.4, ls=":")
ax.text(6.55, 4.9, "sees c, knows Enc and Dec;\nnot k", fontsize=8.8, color=RED, va="center")
ax.annotate("", xy=(10.8, 7.05), xytext=(1.2, 7.05),
            arrowprops=dict(arrowstyle="<|-|>", color=PURPLE, lw=1.4, ls="--",
                            connectionstyle="arc3,rad=-0.18", shrinkA=4, shrinkB=4))
ax.text(6.0, 8.0, "key k delivered in advance (courier, factory, a meeting)", fontsize=8.8,
        color=PURPLE, ha="center")

ax.plot([0.2, 11.8], [3.75, 3.75], color="#d1d5db", lw=1)
ax.text(0.2, 3.35, "(b) public key: Bob publishes pk, keeps sk", fontsize=11, weight="bold", color=INK)
y = 1.9
node(ax, 1.0, y, "Alice", color=BLUE, r=0.55, fontsize=9.5)
box(ax, 3.4, y, "Enc$_{pk}(m)$", w=1.8, h=0.8, color=BLUE)
box(ax, 8.6, y, "Dec$_{sk}(c)$", w=1.8, h=0.8, color=GREEN)
node(ax, 11.0, y, "Bob", color=GREEN, r=0.55, fontsize=9.5)
arrow(ax, 1.0, y, 3.4, y, text="m", shrink=30)
arrow(ax, 3.4, y, 8.6, y, text="ciphertext c", shrink=45, color=INK)
arrow(ax, 8.6, y, 11.0, y, text="m", shrink=40)
box(ax, 6.0, 0.45, "public directory: pk", w=2.6, h=0.6, color=ORANGE, fc="#fff7ed")
arrow(ax, 10.6, 1.45, 7.3, 0.55, color=ORANGE, shrink=6, curve=-0.2)
arrow(ax, 4.7, 0.55, 3.4, 1.45, color=ORANGE, shrink=6, curve=-0.2)
node(ax, 6.0, 2.95, "Eve", color=RED, r=0.45, fontsize=9.5)
ax.plot([6.0, 6.0], [1.9, 2.5], color=RED, lw=1.4, ls=":")
ax.text(6.55, 2.95, "sees c and pk", fontsize=8.8, color=RED, va="center")
ax.text(11.0, 0.95, "sk never\nleaves Bob", fontsize=8.8, color=GREEN, ha="center")
save(fig, "ch77-alice-bob-eve")

# ---------------------------------------------------------------- figure 2: toy RSA flow
p, q, e = 101, 113, 3
n, phi = p * q, (p - 1) * (q - 1)
d = pow(e, -1, phi)
m = 4242
c = pow(m, e, n)
assert pow(c, d, n) == m

fig, ax = diagram(figsize=(9.6, 4.6))
ax.set_xlim(0, 12); ax.set_ylim(0, 5.8)
box(ax, 2.0, 4.5, f"KEY GENERATION (Bob)\nprimes p = {p}, q = {q}\nn = pq = {n}\n"
    f"φ(n) = (p−1)(q−1) = {phi}", w=3.6, h=1.9, color=PURPLE, fontsize=9.2)
box(ax, 6.4, 5.0, f"public key (n, e) = ({n}, {e})", w=3.4, h=0.7, color=ORANGE, fc="#fff7ed")
box(ax, 7.0, 3.9, f"private key d = e⁻¹ mod φ(n) = {d}\n(check: {e}·{d} = {e*d} ≡ 1)",
    w=3.9, h=0.95, color=GREEN, fc="#f0fdf4", fontsize=9)
arrow(ax, 3.8, 4.9, 4.7, 5.0, shrink=2)
arrow(ax, 3.8, 4.1, 5.05, 3.95, shrink=2)
box(ax, 2.0, 1.4, f"ENCRYPT (Alice)\nmessage m = {m}\nc = mᵉ mod n = {c}", w=3.4, h=1.4, color=BLUE, fontsize=9.3)
box(ax, 9.8, 1.4, f"DECRYPT (Bob)\nc = {c}\ncᵈ mod n = {pow(c, d, n)} = m", w=3.4, h=1.4, color=GREEN, fontsize=9.3)
arrow(ax, 3.7, 1.4, 8.1, 1.4, text=f"c = {c} over the open network", shrink=4, color=INK)
arrow(ax, 4.9, 4.65, 2.8, 2.1, shrink=4, color=ORANGE, curve=0.1)
arrow(ax, 8.95, 3.9, 9.8, 2.1, shrink=4, color=GREEN, curve=-0.15)
ax.text(10.1, 4.55, "Eve must factor\nn = 11413 to get φ(n)\n(trivial here; infeasible\nfor 2048-bit n)",
        fontsize=8.6, color=RED, ha="left", va="center")
save(fig, "ch77-toy-rsa")

# ---------------------------------------------------------------- figure 3: ECB vs CTR with a toy block cipher
H, W = 48, 64
yy, xx = np.mgrid[0:H, 0:W]
img = np.full((H, W), 2, dtype=int)                       # background level 2
img[((xx - 20) ** 2 + (yy - 24) ** 2) < 14 ** 2] = 9      # a disc
img[(xx > 38) & (xx < 58) & (yy > 8) & (yy < 40)] = 14    # a square
img[(xx > 42) & (xx < 54) & (yy > 18) & (yy < 30)] = 5    # hole in the square
rng = np.random.default_rng(77)
perm = rng.permutation(1 << 16)                           # the "key": a random 16-bit permutation


def blocks(a):
    """four 4-bit pixels -> one 16-bit block"""
    a = a.reshape(H, W // 4, 4)
    return (a[..., 0] << 12) | (a[..., 1] << 8) | (a[..., 2] << 4) | a[..., 3]


def pixels(b):
    return np.stack([(b >> 12) & 15, (b >> 8) & 15, (b >> 4) & 15, b & 15], axis=-1).reshape(H, W)


B = blocks(img)
ecb = pixels(perm[B])
ctr_idx = (np.arange(B.size).reshape(B.shape) * 7919 + 12345) % (1 << 16)
ctr = pixels(B ^ perm[ctr_idx])

fig, axs = new_figure(1, 3, figsize=(10.5, 3.3))
for a, im, t in zip(axs, [img, ecb, ctr], ["plaintext image", "ECB: each block encrypted alone",
                                            "CTR: block ⊕ Enc(counter)"]):
    a.imshow(im, cmap="viridis", vmin=0, vmax=15, interpolation="nearest")
    a.set_title(t, fontsize=10.5); a.set_xticks([]); a.set_yticks([]); a.grid(False)
save(fig, "ch77-ecb-vs-ctr")
