import sys, pathlib; sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
from style import plt, np, PALETTE, INK, MUTED, new_figure, diagram, node, box, arrow, save
import scipy.sparse as sp

rng = np.random.default_rng(87)
N, KMEAN = 20000, 6


def er_edges(n, kmean):
    """Erdos-Renyi G(n, p) with p = kmean/(n-1), sampled as a random edge list."""
    m = rng.binomial(n * (n - 1) // 2, kmean / (n - 1))
    u, v = rng.integers(0, n, m), rng.integers(0, n, m)
    keep = u != v
    return u[keep], v[keep]


def ba_edges(n, m):
    """Barabasi-Albert: each new node attaches m edges, preferentially by degree."""
    targets = list(range(m))
    repeated = []                      # every node appears once per unit of degree
    us, vs = [], []
    for new in range(m, n):
        for t in set(targets):
            us.append(new); vs.append(t)
        repeated.extend(targets)
        repeated.extend([new] * m)
        chosen = set()
        while len(chosen) < m:
            chosen.add(repeated[rng.integers(len(repeated))])
        targets = list(chosen)
    return np.array(us), np.array(vs)


def adjacency(n, u, v):
    A = sp.coo_matrix((np.ones(len(u)), (u, v)), shape=(n, n)).tocsr()
    A = ((A + A.T) > 0).astype(float)
    return A


A_er = adjacency(N, *er_edges(N, KMEAN))
A_ba = adjacency(N, *ba_edges(N, KMEAN // 2))
deg_er = np.asarray(A_er.sum(1)).ravel()
deg_ba = np.asarray(A_ba.sum(1)).ravel()

# --- Figure 1: degree distributions on log-log axes -----------------------------------
fig, ax = new_figure(figsize=(7.0, 4.4))
bins = np.unique(np.round(np.logspace(0, np.log10(deg_ba.max() + 1), 30)))
for deg, col, lab in [(deg_er, PALETTE[0], "Erdős–Rényi (random)"),
                      (deg_ba, PALETTE[1], "Barabási–Albert (preferential attachment)")]:
    h, e = np.histogram(deg, bins=bins)
    w = np.diff(e)
    c = np.sqrt(e[:-1] * e[1:])
    ok = h > 0
    ax.loglog(c[ok], h[ok] / w[ok] / len(deg), "o", color=col, ms=5, label=lab)
from scipy.stats import poisson
kk = np.arange(1, 20)
ax.loglog(kk, poisson.pmf(kk, deg_er.mean()), color=PALETTE[0], lw=1.2, alpha=0.7)
kk = np.logspace(np.log10(3), np.log10(deg_ba.max()), 50)
ax.loglog(kk, 2 * 3 ** 2 / kk ** 3, "--", color=PALETTE[1], lw=1.4)
ax.text(60, 2 * 9 / 60 ** 3 * 2.2, r"$P(k)=2m^2/k^3$", color=PALETTE[1], fontsize=10)
ax.text(1.05, 2.5e-3, "Poisson,\nmean 6", color=PALETTE[0], fontsize=10)
ax.set_xlabel("degree $k$")
ax.set_ylabel("fraction of nodes $P(k)$")
ax.set_ylim(1e-7, 1)
ax.legend(loc="lower left")
ax.set_title(f"Same size, same mean degree, different worlds (N = {N})")
save(fig, "ch87-degree-distributions")

# --- Figure 2: SIR on ER vs BA --------------------------------------------------------
beta, gamma, T = 0.05, 0.2, 250


def sir(A, seeds=10):
    n = A.shape[0]
    S = np.ones(n, bool); I = np.zeros(n, bool)
    idx = rng.choice(n, seeds, replace=False)
    S[idx] = False; I[idx] = True
    curve = []
    for _ in range(T):
        k_inf = A @ I.astype(float)
        newI = S & (rng.random(n) < 1 - (1 - beta) ** k_inf)
        rec = I & (rng.random(n) < gamma)
        S &= ~newI; I = (I | newI) & ~rec
        curve.append(I.mean())
    return np.array(curve), 1 - S.mean()


fig, ax = new_figure(figsize=(7.0, 4.0))
t = np.arange(1, T + 1)
for A, col, lab in [(A_er, PALETTE[0], "Erdős–Rényi"), (A_ba, PALETTE[1], "Barabási–Albert")]:
    c, final = sir(A)
    ax.plot(t, 100 * c, color=col, label=f"{lab}: {100 * final:.0f}% ever infected")
ax.set_xlabel("time step")
ax.set_ylabel("currently infected (%)")
ax.set_title(r"SIR epidemic, transmission $\beta=0.05$ per contact, recovery $\gamma=0.2$")
ax.legend(loc="upper right")
save(fig, "ch87-sir-networks")

# --- Figure 3: Kuramoto order parameter vs coupling -----------------------------------
n, g = 2000, 0.5
omega = g * np.tan(np.pi * (np.arange(1, n + 1) - 0.5) / n - np.pi / 2)   # Lorentzian quantiles
Ks = np.linspace(0, 3, 25)
rs = []
for K in Ks:
    th = rng.uniform(0, 2 * np.pi, n)
    dt, rr = 0.05, []
    for step in range(4000):
        z = np.exp(1j * th).mean()
        th = th + dt * (omega + K * np.abs(z) * np.sin(np.angle(z) - th))
        if step >= 2000:
            rr.append(np.abs(z))
    rs.append(np.mean(rr))
fig, ax = new_figure(figsize=(7.0, 4.0))
Kt = np.linspace(1, 3, 200)
ax.plot(Kt, np.sqrt(1 - 2 * g / Kt), color=INK, lw=1.5, label=r"theory $r=\sqrt{1-K_c/K}$")
ax.plot([0, 1], [0, 0], color=INK, lw=1.5)
ax.plot(Ks, rs, "o", color=PALETTE[3], ms=6, label=f"simulation, N = {n}")
ax.axvline(2 * g, color=MUTED, ls=":", lw=1.2)
ax.text(2 * g + 0.04, 0.85, r"$K_c=2\gamma=1$", color=MUTED)
ax.set_xlabel("coupling strength $K$")
ax.set_ylabel("order parameter $r$")
ax.set_ylim(-0.03, 1)
ax.set_title("Kuramoto oscillators: a sudden onset of synchrony")
ax.legend(loc="lower right")
save(fig, "ch87-kuramoto")
print("ER mean deg", deg_er.mean(), "max", deg_er.max(), "BA mean", deg_ba.mean(), "max", deg_ba.max())
print("Kuramoto", list(zip(Ks.round(2), np.round(rs, 3))))
