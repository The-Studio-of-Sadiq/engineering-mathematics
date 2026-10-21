"""Chapter 40 figures: the Prisoner's Dilemma matrix; a repeated-game round-robin tournament."""
import sys, pathlib; sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
from style import plt, np, PALETTE, INK, MUTED, new_figure, diagram, node, box, arrow, save
from matplotlib.patches import Rectangle, Ellipse

BLUE, RED, GREEN, GREY = PALETTE[0], PALETTE[1], PALETTE[2], PALETTE[6]

# ---------- Figure 1: the matrix ----------
fig, ax = plt.subplots(figsize=(8.4, 4.4), constrained_layout=True)
ax.set_xlim(-1.6, 9.2); ax.set_ylim(-3.9, 1.5); ax.set_aspect("equal"); ax.axis("off")
cell, x0, y0 = 1.8, 0.4, 0.0
A = np.array([[3, 0], [5, 1]]); B = A.T
rows, cols = ["Cooperate", "Defect"], ["Cooperate", "Defect"]
lab = [["R , R", "S , T"], ["T , S", "P , P"]]
ax.text(x0 + cell, y0 + 1.0, "Player 2", ha="center", fontsize=10.5, color=RED, weight="bold")
ax.text(x0 - 1.45, y0 - cell, "Player 1", ha="center", va="center", rotation=90, fontsize=10.5, color=BLUE, weight="bold")
for j, c in enumerate(cols):
    ax.text(x0 + (j + 0.5) * cell, y0 + 0.35, c, ha="center", va="center", fontsize=10, color=RED)
for i, r in enumerate(rows):
    ax.text(x0 - 0.55, y0 - (i + 0.5) * cell, r, ha="center", va="center", fontsize=10, color=BLUE, rotation=90)
for i in range(2):
    for j in range(2):
        x, y = x0 + j * cell, y0 - (i + 1) * cell
        fc = "#fee2e2" if (i, j) == (1, 1) else ("#dcfce7" if (i, j) == (0, 0) else "white")
        ax.add_patch(Rectangle((x, y), cell, cell, fc=fc, ec=INK, lw=1.2))
        xa, xb, yc = x + 0.3 * cell, x + 0.7 * cell, y + 0.58 * cell
        ax.text(xa, yc, f"{A[i, j]}", ha="center", va="center", fontsize=13, color=BLUE)
        ax.text(x + 0.5 * cell, yc, ",", ha="center", va="center", fontsize=13, color=INK)
        ax.text(xb, yc, f"{B[i, j]}", ha="center", va="center", fontsize=13, color=RED)
        ax.text(x + 0.5 * cell, y + 0.22 * cell, lab[i][j], ha="center", va="center", fontsize=9, color=MUTED)
        if A[i, j] == A[:, j].max():
            ax.add_patch(Ellipse((xa, yc), 0.3 * cell, 0.3 * cell, fc="none", ec=BLUE, lw=1.4))
        if B[i, j] == B[i, :].max():
            ax.add_patch(Ellipse((xb, yc), 0.3 * cell, 0.3 * cell, fc="none", ec=RED, lw=1.4))
tx = 5.0
ax.text(tx, 0.9, "Why it is a trap", fontsize=11.5, weight="bold", color=INK, va="top")
notes = [("T > R", "temptation beats reward", INK),
         ("P > S", "punishment beats being the sucker", INK),
         ("R > P", "mutual cooperation beats mutual defection", INK),
         ("", "", INK),
         ("green", "what both would choose together", GREEN),
         ("red", "the only Nash equilibrium", RED),
         ("circles", "best responses: Defect, every time", MUTED)]
for k, (a, b, c) in enumerate(notes):
    if a:
        ax.text(tx, 0.25 - 0.52 * k, a, fontsize=10, color=c, weight="bold", va="center")
        ax.text(tx + 1.25, 0.25 - 0.52 * k, b, fontsize=9.5, color=INK, va="center")
save(fig, "ch40-pd-matrix")

# ---------- Figure 2: tournament ----------
T, R, P, S = 5, 3, 1, 0
PAY = {('C', 'C'): (R, R), ('C', 'D'): (S, T), ('D', 'C'): (T, S), ('D', 'D'): (P, P)}
rng = np.random.default_rng(7)
STRATS = {
    'ALL-C': lambda me, op: 'C',
    'ALL-D': lambda me, op: 'D',
    'TFT': lambda me, op: op[-1] if op else 'C',
    'GRIM': lambda me, op: 'D' if 'D' in op else 'C',
    'TF2T': lambda me, op: 'D' if op[-2:] == ['D', 'D'] else 'C',
    'WSLS': lambda me, op: 'C' if not me else (me[-1] if op[-1] == 'C' else ('D' if me[-1] == 'C' else 'C')),
    'RANDOM': lambda me, op: 'C' if rng.random() < 0.5 else 'D',
}


def play(f, g, rounds, noise):
    a, b, sa, sb = [], [], 0, 0
    for _ in range(rounds):
        x, y = f(a, b), g(b, a)
        if rng.random() < noise: x = 'D' if x == 'C' else 'C'
        if rng.random() < noise: y = 'D' if y == 'C' else 'C'
        pa, pb = PAY[(x, y)]
        a.append(x); b.append(y); sa += pa; sb += pb
    return sa / rounds, sb / rounds


def tournament(noise, reps, rounds=200):
    names = list(STRATS)
    tot = {n: 0.0 for n in names}
    for _ in range(reps):
        score = {n: [] for n in names}
        for i, n in enumerate(names):
            for m in names[i:]:
                s1, s2 = play(STRATS[n], STRATS[m], rounds, noise)
                score[n].append(s1)
                if m != n: score[m].append(s2)
        for n in names: tot[n] += np.mean(score[n]) / reps
    return tot


clean, noisy = tournament(0.0, 20), tournament(0.05, 20)
order = sorted(clean, key=lambda n: -clean[n])
print({n: round(clean[n], 3) for n in order}); print({n: round(noisy[n], 3) for n in order})
fig, ax = new_figure(figsize=(7.6, 4.2))
xs = np.arange(len(order)); w = 0.38
b1 = ax.bar(xs - w / 2, [clean[n] for n in order], w, color=BLUE, label="no noise")
b2 = ax.bar(xs + w / 2, [noisy[n] for n in order], w, color=PALETTE[4], label="5% of moves flipped")
for bars in (b1, b2):
    for r in bars:
        ax.text(r.get_x() + r.get_width() / 2, r.get_height() + 0.03, f"{r.get_height():.2f}", ha="center",
                va="bottom", fontsize=7.5, color=INK)
ax.axhline(R, color=GREEN, ls="--", lw=1.2, label="R = 3: two cooperators")
ax.axhline(P, color=RED, ls="--", lw=1.2, label="P = 1: two defectors")
ax.set_xticks(xs); ax.set_xticklabels(order)
ax.set_ylim(0, 3.9); ax.set_ylabel("average payoff per round")
ax.set_title("Round-robin tournament, 200 rounds per match")
ax.legend(loc="upper center", ncol=4, bbox_to_anchor=(0.5, 1.0), fontsize=8.5)
ax.grid(axis="x", visible=False)
save(fig, "ch40-tournament")
