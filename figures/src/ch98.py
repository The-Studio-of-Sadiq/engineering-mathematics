import sys, pathlib; sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
from style import plt, np, PALETTE, INK, MUTED, new_figure, diagram, node, box, arrow, save

# --- Figure 1: the conjecture lifecycle -------------------------------------------------
fig, ax = diagram(figsize=(8.6, 5.0))
ax.set_xlim(0, 11); ax.set_ylim(0, 6.4)
B = dict(w=2.5, h=0.95, fontsize=9.5)
box(ax, 1.6, 5.2, "Observation\n(data, pattern, analogy)", color=PALETTE[6], **B)
box(ax, 5.5, 5.2, "Conjecture\n(precise statement)", color=PALETTE[0], **B)
box(ax, 9.4, 5.2, "Evidence\ncomputation, heuristics,\nspecial cases", color=PALETTE[5], w=2.6, h=1.1, fontsize=9.2)
box(ax, 9.4, 2.6, "Counterexample\nsearch", color=PALETTE[4], **B)
box(ax, 5.5, 1.0, "Proof\n(theorem)", color=PALETTE[2], **B)
box(ax, 9.4, 0.55, "Refutation\n(counterexample)", color=PALETTE[1], w=2.6, h=0.85, fontsize=9.2)
box(ax, 1.6, 2.6, "New questions\n(sharper, stronger,\ngeneralised)", color=PALETTE[3], w=2.6, h=1.1, fontsize=9.2)
box(ax, 5.5, 3.1, "Proof strategy\nreduce, induct, symmetrise,\nlet a machine check cases", color=PALETTE[0], w=3.0, h=1.1,
    fontsize=9.0, fc="#eef3ff")
def ar(x0, y0, x1, y1, col=MUTED, curve=0.0):
    ax.annotate("", xy=(x1, y1), xytext=(x0, y0), zorder=5,
                arrowprops=dict(arrowstyle="-|>", color=col, lw=1.6, shrinkA=1, shrinkB=1,
                                connectionstyle=f"arc3,rad={curve}"))
ar(2.85, 5.2, 4.25, 5.2); ax.text(3.55, 5.32, "guess", fontsize=9, color=MUTED, ha="center")
ar(6.75, 5.2, 8.1, 5.2); ax.text(7.42, 5.32, "test", fontsize=9, color=MUTED, ha="center")
ar(9.4, 4.65, 9.4, 3.08); ax.text(9.55, 3.9, "try to break it", fontsize=9, color=MUTED, va="center")
ar(9.4, 2.12, 9.4, 0.98); ax.text(9.55, 1.55, "found one", fontsize=9, color=PALETTE[1], va="center")
ar(8.15, 2.7, 7.0, 2.95); ax.text(7.6, 2.4, "survives", fontsize=9, color=MUTED, ha="center")
ar(5.5, 2.55, 5.5, 1.48)
ar(5.5, 3.65, 5.5, 4.72, col=PALETTE[0]); ax.text(5.65, 4.2, "modify the\nstatement", fontsize=8.5, color=PALETTE[0], va="center")
ar(4.25, 1.15, 2.6, 2.05)
ar(8.1, 0.4, 1.6, 2.05, col=PALETTE[1], curve=-0.18)
ax.text(4.0, 0.05, "a refutation is new knowledge too", fontsize=8.8, color=PALETTE[1], ha="center")
ar(1.6, 3.15, 1.6, 4.72); ax.text(1.45, 3.95, "repeat", fontsize=9, color=MUTED, ha="right", va="center")
save(fig, "ch98-conjecture-lifecycle")

# --- Figure 2: two counterexamples that came late ---------------------------------------
N = 2_000_000
s = np.ones(N + 1, bool); s[:2] = False
for i in range(2, int(N ** 0.5) + 1):
    if s[i]:
        s[i * i::i] = False
p = np.nonzero(s)[0]
r = np.zeros(N + 1, int); r[p[p % 4 == 3]] += 1; r[p[p % 4 == 1]] -= 1
d = np.cumsum(r)                                     # pi(x;4,3) - pi(x;4,1)
mu = np.ones(N + 1, int); mu[0] = 0
for q in p:
    mu[q::q] *= -1
    if q * q <= N:
        mu[q * q::q * q] = 0
M = np.cumsum(mu)
x = np.arange(N + 1)

fig, (a1, a2) = new_figure(1, 2, figsize=(10.5, 4.0))
idx = np.unique(np.geomspace(3, N, 6000).astype(int))
a1.plot(x[idx], d[idx], color=PALETTE[0], lw=1.2)
a1.axhline(0, color=INK, lw=0.8)
for x0 in (26861, 616841):
    a1.axvline(x0, color=PALETTE[1], lw=1, ls="--")
a1.text(5, 145, "Team 4k+3 leads almost everywhere.\nDashed: 4k+1 first leads at x = 26,861,\nthen again at x = 616,841.",
        va="top", fontsize=9, color=PALETTE[1])
a1.set_xscale("log"); a1.set_xlim(3, N)
a1.set_xlabel("x"); a1.set_ylabel(r"$\pi(x;4,3)-\pi(x;4,1)$")
a1.set_title("The prime race mod 4")

xs = np.unique(np.geomspace(1, N, 8000).astype(int))
a2.plot(xs, M[xs] / np.sqrt(xs), color=PALETTE[3], lw=1.0, label=r"$M(x)/\sqrt{x}$")
a2.axhline(1, color=PALETTE[1], ls="--", lw=1.2)
a2.text(30, 1.06, "Mertens' guess: |M(x)| < sqrt(x) forever", color=PALETTE[1], fontsize=9)
a2.axhline(-1, color=PALETTE[1], ls="--", lw=1.2)
a2.set_xscale("log"); a2.set_xlim(1, N); a2.set_ylim(-1.25, 1.25)
a2.set_xlabel("x"); a2.set_ylabel(r"$M(x)/\sqrt{x}$")
a2.set_title("Mertens' conjecture: false, yet never seen to fail")
save(fig, "ch98-late-counterexamples")

# --- Figure 3: map of open problems referenced across the book --------------------------
tracks = [
    ("A  Algorithms", PALETTE[0], ["Is ω = 2? (Ch. 26)", "R(5,5): 43–46 (Ch. 28)", "P vs NP (Ch. 30)",
                                    "Collatz (Ch. 29, 30)"]),
    ("B  Probability", PALETTE[5], ["Ramsey growth base:\n√2 or near 4? (Ch. 32)"]),
    ("C  Game theory", PALETTE[4], ["Nash equilibria fast?\nPPAD vs P (Ch. 39)"]),
    ("D  Optimization", PALETTE[2], ["Why non-convex training\nworks (Ch. 47, 53)", "Why real integer\nprograms are easy (Ch. 49)"]),
    ("E  Abstract", PALETTE[3], ["Riemann hypothesis\n(Ch. 57)", "Twin primes (Ch. 57)"]),
    ("F  Logic", PALETTE[6], ["Goldbach (Ch. 63)", "Why SAT solvers\nwork (Ch. 67)"]),
    ("G  Graphics", PALETTE[1], ["MLC: Mandelbrot\nlocally connected? (Ch. 73)"]),
    ("H  Crypto", PALETTE[0], ["Lattices vs quantum\ncomputers (Ch. 77)", "One-way functions\nexist? (Ch. 77)"]),
    ("I  Dynamics", PALETTE[5], ["Is rule 30 random?\n(Ch. 86)"]),
]
fig, ax = plt.subplots(figsize=(10.5, 5.6), constrained_layout=True)
ax.set_xlim(0, 10.5); ax.set_ylim(-0.4, 9.6); ax.axis("off")
for i, (name, col, probs) in enumerate(tracks):
    y = 9.0 - i * 1.05
    box(ax, 1.05, y, name, w=1.9, h=0.7, color=col, fc="white", fontsize=9.5)
    for k, pr in enumerate(probs):
        xk = 3.15 + k * 2.1
        box(ax, xk, y, pr, w=2.0, h=0.86, color=col, fc="#f8fafc", fontsize=7.9, lw=1.2)
        ax.plot([2.0 if k == 0 else xk - 2.1 + 1.0, xk - 1.0], [y, y], color=col, lw=1.0, zorder=1)
ax.text(10.4, -0.25, "Solved since posed: four-colour theorem (1976), Kepler (1998), R(4,4) = 18 (1955)",
        ha="right", fontsize=8.8, color=MUTED)
ax.set_title("Open problems this book has pointed at, by track", loc="left")
save(fig, "ch98-open-problems-map")
