"""Figures for Chapter 79 — Error-Correcting Codes: Recovering the Lost."""
import sys, pathlib; sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
from style import plt, np, PALETTE, INK, MUTED, new_figure, diagram, node, box, arrow, save
from matplotlib.patches import Circle

BLUE, RED, GREEN, PURPLE, ORANGE, TEAL, GREY = PALETTE

# ---------------------------------------------------------------- figure 1: codewords on the cube
# oblique projection of the unit cube: bit order (b1 b2 b3) -> x, y, depth
def proj(v):
    x, y, z = v
    return np.array([2.2 * x + 1.1 * z, 2.2 * y + 0.75 * z])


verts = [(a, b, c) for a in (0, 1) for b in (0, 1) for c in (0, 1)]
edges = [(u, v) for u in verts for v in verts if u < v and sum(abs(i - j) for i, j in zip(u, v)) == 1]


def draw_cube(ax, colour_of, ring_of, title):
    for u, v in edges:
        pu, pv = proj(u), proj(v)
        ax.plot([pu[0], pv[0]], [pu[1], pv[1]], color="#cbd5e1", lw=1.6, zorder=1)
    for v in verts:
        pv = proj(v)
        c = colour_of(v)
        ring = ring_of(v)
        if ring:
            ax.add_patch(Circle(pv, 0.42, fc=c, ec=c, lw=2.4, zorder=3))
            ax.text(*pv, "".join(map(str, v)), ha="center", va="center", fontsize=10, color="white",
                    weight="bold", zorder=4)
        else:
            ax.add_patch(Circle(pv, 0.36, fc="white", ec=c, lw=2.0, zorder=3))
            ax.text(*pv, "".join(map(str, v)), ha="center", va="center", fontsize=9.5, color=INK, zorder=4)
    ax.set_xlim(-0.8, 4.2); ax.set_ylim(-0.8, 3.8); ax.set_aspect("equal"); ax.axis("off")
    ax.set_title(title, fontsize=10.5)


fig, axs = new_figure(1, 2, figsize=(10.0, 4.6))
w = lambda v: sum(v)
draw_cube(axs[0], lambda v: BLUE if w(v) <= 1 else RED, lambda v: v in [(0, 0, 0), (1, 1, 1)],
          "repetition code {000, 111}: distance 3\nblue sphere decodes to 000, red to 111")
draw_cube(axs[1], lambda v: GREEN if w(v) % 2 == 0 else GREY, lambda v: w(v) % 2 == 0,
          "even-parity code {000, 011, 101, 110}: distance 2\nevery error is seen, none can be located")
save(fig, "ch79-hamming-cube")

# ---------------------------------------------------------------- figure 2: Hamming(7,4) Venn diagram
centres = {"A": (0.0, 0.62), "B": (0.72, 0.62), "C": (0.36, 0.0)}
R = 0.78
# region label positions: positions 1..7 = p1 p2 d1 p3 d2 d3 d4
spots = {1: (-0.42, 0.95), 2: (1.14, 0.95), 3: (0.36, 0.92), 4: (0.36, -0.42),
         5: (-0.08, 0.22), 6: (0.80, 0.22), 7: (0.36, 0.40)}
names = {1: "p₁", 2: "p₂", 3: "d₁", 4: "p₃", 5: "d₂", 6: "d₃", 7: "d₄"}
code = [0, 1, 1, 0, 0, 1, 1]               # data 1011 -> codeword p1 p2 d1 p3 d2 d3 d4
recv = code.copy(); recv[5] ^= 1            # flip position 6 (d3)
members = {"A": [1, 3, 5, 7], "B": [2, 3, 6, 7], "C": [4, 5, 6, 7]}
cols = {"A": BLUE, "B": ORANGE, "C": GREEN}

fig, axs = new_figure(1, 2, figsize=(10.4, 4.8))
for ax, word, title in [(axs[0], code, "encode data 1011: each circle gets even parity"),
                        (axs[1], recv, "one bit flipped: circles B and C fail,\nA passes, so the error is in B ∩ C, not A")]:
    for k, (cx, cy) in centres.items():
        odd = sum(word[i - 1] for i in members[k]) % 2
        ax.add_patch(Circle((cx, cy), R, fc=cols[k], alpha=0.07, ec=RED if odd else cols[k],
                            lw=3.0 if odd else 1.8, ls="-" if not odd else "--"))
        lx, ly = {"A": (-0.75, 1.3), "B": (1.47, 1.3), "C": (0.36, -0.98)}[k]
        ax.text(lx, ly, f"{k}: {'odd!' if odd else 'even'}", color=RED if odd else cols[k], fontsize=10.5,
                ha="center", weight="bold")
    for pos, (x, y) in spots.items():
        bad = word is recv and pos == 6
        ax.text(x, y + 0.09, str(word[pos - 1]), ha="center", va="center", fontsize=17,
                color=RED if bad else INK, weight="bold")
        ax.text(x, y - 0.13, names[pos], ha="center", va="center", fontsize=9, color=MUTED)
    ax.set_xlim(-1.05, 1.8); ax.set_ylim(-1.1, 1.5); ax.set_aspect("equal"); ax.axis("off")
    ax.set_title(title, fontsize=10.5)
save(fig, "ch79-hamming74-venn")

# ---------------------------------------------------------------- figure 3: Reed-Solomon as a polynomial through points
coef = np.array([2.0, 1.0, -0.6, 0.08])     # message (k = 4): f(x) = 2 + x - 0.6x^2 + 0.08x^3
f = lambda x: np.polyval(coef[::-1], x)
xs = np.arange(1, 8)                         # n = 7 evaluation points
ys = f(xs)
erased, corrupt = 3, 6
fig, ax = new_figure(figsize=(8.0, 4.4))
xx = np.linspace(0.5, 7.5, 300)
ax.plot(xx, f(xx), color=BLUE, lw=2.2, label="message polynomial $f(x)$, degree $< k=4$")
ok = [i for i in xs if i not in (erased, corrupt)]
ax.plot(ok, f(np.array(ok)), "o", color=BLUE, ms=9, label="received symbols $f(1),\\dots,f(7)$")
ax.plot([erased], [f(erased)], "o", ms=12, mfc="white", mec=GREY, mew=2, ls="none", label="erased (position known)")
ax.text(erased, f(erased) + 0.45, "?", ha="center", fontsize=15, color=GREY, weight="bold")
yc = f(corrupt) + 1.6
ax.plot([corrupt], [yc], "X", color=RED, ms=12, label="corrupted (position unknown)")
ax.annotate("", xy=(corrupt, f(corrupt) + 0.12), xytext=(corrupt, yc - 0.15),
            arrowprops=dict(arrowstyle="-|>", color=RED, lw=1.3, ls="--"))
ax.text(corrupt - 0.2, yc + 0.05, "off the curve", color=RED, fontsize=9.5, ha="right")
ax.set_xticks(xs); ax.set_xlabel("evaluation point $x$"); ax.set_ylabel("symbol value")
ax.legend(loc="upper left", fontsize=9)
ax.set_title("$n=7$ values of a cubic: any 4 correct ones determine it; $d=n-k+1=4$", fontsize=10.5)
ax.set_ylim(0, 8.2)
save(fig, "ch79-reed-solomon")
