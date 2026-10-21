"""Figures for Chapter 78 — Elliptic Curves and Modern Cryptographic Mathematics."""
import sys, pathlib; sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
from style import plt, np, PALETTE, INK, MUTED, new_figure, diagram, node, box, arrow, save

BLUE, RED, GREEN, PURPLE, ORANGE, TEAL, GREY = PALETTE

# ---------------------------------------------------------------- figure 1: chord and tangent on y^2 = x^3 + 17
B = 17
x0 = -(B ** (1 / 3))
xs = np.linspace(x0, 6.2, 1200)
ys = np.sqrt(np.maximum(xs ** 3 + B, 0))
ys[0] = 0.0

fig, axs = new_figure(1, 2, figsize=(10.5, 4.9))
for ax in axs:
    ax.plot(xs, ys, color=BLUE, lw=2.2)
    ax.plot(xs, -ys, color=BLUE, lw=2.2)
    ax.axhline(0, color=MUTED, lw=0.8); ax.axvline(0, color=MUTED, lw=0.8)
    ax.set_xlim(-3.6, 6.2); ax.set_ylim(-11, 11)
    ax.set_xlabel("$x$"); ax.set_ylabel("$y$")

# (a) chord: P=(-2,3), Q=(2,5) -> third point (1/4, 33/8), sum (1/4, -33/8)
ax = axs[0]
P, Q, R = (-2, 3), (2, 5), (0.25, 33 / 8)
lx = np.linspace(-3.6, 6.2, 10)
ax.plot(lx, 0.5 * lx + 4, color=ORANGE, lw=1.6)
ax.plot([0.25, 0.25], [33 / 8, -33 / 8], color=GREY, lw=1.3, ls="--")
for pt, lab, c, off in [(P, "$P=(-2,3)$", INK, (-1.5, 0.9)), (Q, "$Q=(2,5)$", INK, (0.25, -1.4)),
                        (R, "$R$", ORANGE, (-0.75, 0.6)), ((0.25, -33 / 8), "$P+Q=(\\frac{1}{4},-\\frac{33}{8})$", RED, (0.35, -1.5))]:
    ax.plot(*pt, "o", color=c, ms=7, zorder=5)
    ax.text(pt[0] + off[0], pt[1] + off[1], lab, color=c, fontsize=10.5)
ax.set_title("chord: the line through P and Q meets the curve at R;\nreflect R to get P + Q", fontsize=10.5)

# (b) tangent: P=(-1,4), slope 3/8 -> 2P = (137/64, -2651/512)
ax = axs[1]
P = (-1, 4)
x2, y2 = 137 / 64, 2651 / 512
ax.plot(lx, 0.375 * (lx + 1) + 4, color=GREEN, lw=1.6)
ax.plot([x2, x2], [y2, -y2], color=GREY, lw=1.3, ls="--")
for pt, lab, c, off in [(P, "$P=(-1,4)$", INK, (-2.3, 0.7)), ((x2, y2), "$R$", GREEN, (-0.3, 0.8)),
                        ((x2, -y2), "$2P\\approx(2.14,-5.18)$", RED, (0.3, -1.3))]:
    ax.plot(*pt, "o", color=c, ms=7, zorder=5)
    ax.text(pt[0] + off[0], pt[1] + off[1], lab, color=c, fontsize=10.5)
ax.set_title("tangent: to double P, use the tangent line at P;\nit meets the curve once more at R", fontsize=10.5)
save(fig, "ch78-chord-tangent")

# ---------------------------------------------------------------- figure 2: the same law over a finite field
def curve_points(p, a, b):
    return [(x, y) for x in range(p) for y in range(p) if (y * y - (x ** 3 + a * x + b)) % p == 0]


def padd(P, Q, p, a):
    if P is None: return Q
    if Q is None: return P
    (x1, y1), (x2, y2) = P, Q
    if x1 == x2 and (y1 + y2) % p == 0: return None
    if P == Q: lam = (3 * x1 * x1 + a) * pow(2 * y1, -1, p) % p
    else: lam = (y2 - y1) * pow(x2 - x1, -1, p) % p
    x3 = (lam * lam - x1 - x2) % p
    return (x3, (lam * (x1 - x3) - y1) % p)


fig, axs = new_figure(1, 2, figsize=(10.8, 5.0))
p, a, b = 17, 2, 2
pts = curve_points(p, a, b)
G = (5, 1)
mult, R = [], None
for k in range(1, 19):
    R = padd(R, G, p, a); mult.append(R)
ax = axs[0]
ax.set_xlim(-0.8, p - 0.2); ax.set_ylim(-0.8, p - 0.2)
ax.set_xticks(range(0, p, 2)); ax.set_yticks(range(0, p, 2)); ax.set_aspect("equal")
ax.axhline(p / 2, color=MUTED, lw=0.9, ls=":")
mx = np.array(mult)
ax.plot(mx[:, 0], mx[:, 1], color=ORANGE, lw=0.9, alpha=0.7, zorder=2)
ax.plot(*np.array(pts).T, "o", color=BLUE, ms=8, zorder=3)
for k, (x, y) in enumerate(mult, start=1):
    ax.text(x + 0.3, y + 0.3, str(k), fontsize=8.5, color=RED if k == 1 else INK, zorder=4)
ax.set_xlabel("$x$ (mod 17)"); ax.set_ylabel("$y$ (mod 17)")
ax.set_title("$y^2=x^3+2x+2$ over $\\mathbb{F}_{17}$: 18 points + $\\mathcal{O}$\n"
             "label $k$ marks $kG$ for $G=(5,1)$; $19G=\\mathcal{O}$", fontsize=10.5)

p2 = 263
pts2 = np.array(curve_points(p2, a, b))
ax = axs[1]
ax.plot(pts2[:, 0], pts2[:, 1], "o", color=BLUE, ms=2.6)
ax.axhline(p2 / 2, color=MUTED, lw=0.9, ls=":")
ax.set_xlim(0, p2); ax.set_ylim(0, p2); ax.set_aspect("equal")
ax.set_xlabel("$x$ (mod 263)"); ax.set_ylabel("$y$ (mod 263)")
ax.set_title(f"same equation over $\\mathbb{{F}}_{{263}}$: {len(pts2) + 1} points\n"
             "(Hasse: within $2\\sqrt{263}\\approx32$ of $264$)", fontsize=10.5)
save(fig, "ch78-finite-field-points")
