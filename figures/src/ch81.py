import sys, pathlib; sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
from style import plt, np, PALETTE, INK, MUTED, new_figure, diagram, node, box, arrow, save

# --- Figure 1: lattice structure of linear congruential generators ---------------------
fig, (a1, a2, a3) = new_figure(1, 3, figsize=(11.0, 3.9))

# (a) a small full-period LCG: consecutive pairs fall on a few lines
m, a, c = 256, 37, 1
x = np.arange(m)
a1.plot(x / m, ((a * x + c) % m) / m, "o", ms=2.6, color=PALETTE[1])
a1.set_title("LCG $x\\mapsto 37x+1$ mod 256", fontsize=10.5)
a1.set_xlabel("$u_k$"); a1.set_ylabel("$u_{k+1}$")
a1.set_xlim(0, 1); a1.set_ylim(0, 1); a1.set_aspect("equal")

# (b) RANDU triples seen edge-on: 15 planes
M, A = 2 ** 31, 65539
xs = np.empty(6000, dtype=np.int64); v = 1
for i in range(xs.size):
    v = (A * v) % M; xs[i] = v
u = xs / M
P = np.stack([u[:-2], u[1:-1], u[2:]], axis=1)
nrm = np.array([9.0, -6.0, 1.0]); nrm /= np.linalg.norm(nrm)
e1 = np.cross(nrm, [0.0, 0.0, 1.0]); e1 /= np.linalg.norm(e1)
a2.plot(P @ nrm, P @ e1, "o", ms=1.4, color=PALETTE[3], alpha=0.6)
a2.set_title("RANDU triples, viewed edge-on", fontsize=10.5)
a2.set_xlabel("position along $(9,-6,1)$")
a2.set_ylabel("position within a plane")
a2.text(0.5, 1.02, "every point lies on one of 15 planes", transform=a2.transAxes,
        ha="center", va="bottom", fontsize=8.8, color=MUTED)
a2.set_title("RANDU triples, viewed edge-on\n", fontsize=10.5)

# (c) a modern generator
rng = np.random.default_rng(2024)
w = rng.random(m + 1)
a3.plot(w[:-1], w[1:], "o", ms=2.6, color=PALETTE[2])
a3.set_title("256 pairs from PCG64", fontsize=10.5)
a3.set_xlabel("$u_k$"); a3.set_ylabel("$u_{k+1}$")
a3.set_xlim(0, 1); a3.set_ylim(0, 1); a3.set_aspect("equal")
save(fig, "ch81-lcg-lattice")

# --- Figure 2: the distinguishing game --------------------------------------------------
fig, ax = diagram(figsize=(9.0, 4.4))
ax.set_xlim(-0.2, 11); ax.set_ylim(0.2, 6.0)
box(ax, 2.7, 3.2, "Challenger\nflips a secret coin $b$", w=3.0, h=1.1, color=PALETTE[0], fontsize=10)
box(ax, 1.4, 5.2, "$b=0$: truly random\n$\\ell$-bit string", w=2.2, h=0.95, color=PALETTE[2], fontsize=9)
box(ax, 4.0, 5.2, "$b=1$: $G(\\text{seed})$,\nseed of $s<\\ell$ bits", w=2.3, h=0.95, color=PALETTE[1], fontsize=9)
arrow(ax, 1.4, 4.72, 2.2, 3.8, shrink=4)
arrow(ax, 4.0, 4.72, 3.2, 3.8, shrink=4)
box(ax, 8.4, 3.2, "Distinguisher $D$\n(time $t$ budget)", w=3.0, h=1.1, color=PALETTE[3], fontsize=10)
arrow(ax, 4.2, 3.2, 6.9, 3.2, text="string $y$", shrink=4, fontsize=10, color=INK)
box(ax, 8.4, 1.2, "outputs a guess $b'$", w=2.6, h=0.75, color=PALETTE[3], fontsize=9.5)
arrow(ax, 8.4, 2.65, 8.4, 1.58, shrink=3)
ax.text(2.9, 1.45, r"advantage $=\left|\,P(b'=1\mid b=1)-P(b'=1\mid b=0)\,\right|$",
        ha="center", fontsize=10, color=INK)
ax.text(2.9, 0.75, r"$G$ is pseudorandom if every efficient $D$ has tiny advantage",
        ha="center", fontsize=9.3, color=MUTED)
save(fig, "ch81-distinguisher-game")
