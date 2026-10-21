import sys, pathlib; sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
from style import plt, np, PALETTE, INK, MUTED, new_figure, diagram, node, box, arrow, save

# --- Figure 1: Monty Hall probability tree (you pick door 1) -------------------
fig, ax = diagram(figsize=(8.6, 5.2))
ax.set_xlim(0, 11.5); ax.set_ylim(0, 6.4)
box(ax, 0.9, 3.2, "You pick\ndoor 1", w=1.5, h=0.9, color=INK)
cars = [(4.0, 5.2, "Car behind 1"), (4.0, 3.2, "Car behind 2"), (4.0, 1.2, "Car behind 3")]
for x, y, t in cars:
    box(ax, x, y, t, w=2.0, h=0.7, color=PALETTE[0])
    arrow(ax, 1.65, 3.2, x - 1.0, y, shrink=4)
    ax.text((1.65 + x - 1.0) / 2 - 0.15, (3.2 + y) / 2 + 0.12, "1/3", color=PALETTE[0], fontsize=10, ha="center")
leaves = [  # (car y, leaf y, label, branch prob, joint, colour)
    (5.2, 5.85, "Host opens 2", "1/2", "1/6", PALETTE[6]),
    (5.2, 4.55, "Host opens 3", "1/2", "1/6", PALETTE[6]),
    (3.2, 3.2, "Host opens 3", "1", "1/3", PALETTE[2]),
    (1.2, 1.2, "Host opens 2", "1", "1/3", PALETTE[2]),
]
for cy, ly, t, bp, joint, col in leaves:
    box(ax, 7.4, ly, t, w=1.9, h=0.55, color=col, fontsize=9)
    arrow(ax, 5.0, cy, 6.45, ly, shrink=3)
    ax.text(5.75, (cy + ly) / 2 + 0.1, bp, fontsize=9.5, color=MUTED, ha="center")
    ax.text(8.55, ly, f"joint {joint}", fontsize=9.5, color=INK, va="center")
    verdict = "stay wins" if col == PALETTE[6] else "switch wins"
    ax.text(10.0, ly, verdict, fontsize=9.5, color=col, va="center", weight="bold")
ax.text(5.75, 0.25, "Given 'host opens 3': stay wins with 1/6, switch wins with 1/3  ->  P(switch wins) = (1/3)/(1/6 + 1/3) = 2/3",
        ha="center", fontsize=9.5, color=INK)
save(fig, "ch31-monty-hall-tree")

# --- Figure 2: running averages, LLN vs Cauchy ----------------------------------
rng = np.random.default_rng(12)
N = 10_000
n = np.arange(1, N + 1)
fig, (a1, a2) = new_figure(1, 2, figsize=(9.2, 3.8))
for k in range(5):
    x = rng.integers(1, 7, N)
    a1.plot(n, np.cumsum(x) / n, lw=1.2, color=PALETTE[k])
a1.axhline(3.5, color=INK, ls="--", lw=1)
band = 3 * np.sqrt(35 / 12) / np.sqrt(n)
a1.fill_between(n, 3.5 - band, 3.5 + band, color=MUTED, alpha=0.15, lw=0, label=r"$\mu \pm 3\sigma/\sqrt{n}$")
a1.set_xscale("log"); a1.set_ylim(1, 6)
a1.set_xlabel("number of rolls $n$"); a1.set_ylabel(r"running average $\bar X_n$")
a1.set_title("Fair die: averages settle at 3.5"); a1.legend(loc="upper right")
for k in range(5):
    x = rng.standard_cauchy(N)
    a2.plot(n, np.cumsum(x) / n, lw=1.2, color=PALETTE[k])
a2.axhline(0, color=INK, ls="--", lw=1)
a2.set_xscale("log"); a2.set_ylim(-8, 8)
a2.set_xlabel("number of samples $n$"); a2.set_ylabel(r"running average $\bar X_n$")
a2.set_title("Cauchy: averages never settle")
save(fig, "ch31-running-averages")
