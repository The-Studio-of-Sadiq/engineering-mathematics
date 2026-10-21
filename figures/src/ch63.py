import sys, pathlib; sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
from style import plt, np, PALETTE, INK, MUTED, new_figure, diagram, node, box, arrow, save
from matplotlib.patches import Rectangle, Polygon
from math import comb

# ---------------------------------------------------------------- 1. proof-strategy flowchart
fig, ax = diagram(figsize=(10.0, 6.4))
ax.set_xlim(0, 13); ax.set_ylim(0, 8.2)

box(ax, 6.5, 7.45, "What is the shape of the claim?", w=4.6, h=0.75, color=INK, fc="#f3f4f6", fontsize=11)

shapes = [
    (1.4, r"$P \Rightarrow Q$", PALETTE[0]),
    (4.0, r"$\forall n\in\mathbb{N}:\ P(n)$", PALETTE[2]),
    (6.6, r"$\exists x:\ P(x)$", PALETTE[3]),
    (9.2, "'no ...', 'not ...',\n'infinitely many'", PALETTE[1]),
    (11.7, r"$\forall x:\ P(x)$" + "\n(but you doubt it)", PALETTE[4]),
]
for x, t, c in shapes:
    box(ax, x, 5.6, t, w=2.3, h=0.95, color=c, fontsize=10)
    arrow(ax, 6.5, 7.05, x, 6.1, shrink=2, lw=1.3)

methods = [
    (1.4, "Direct proof\nassume $P$, walk to $Q$", PALETTE[0]),
    (4.0, "Induction\nbase case + step", PALETTE[2]),
    (6.6, "Construction\nbuild the $x$", PALETTE[3]),
    (9.2, "Contradiction\nassume the opposite", PALETTE[1]),
    (11.7, "Counterexample\nhunt for one $x$", PALETTE[4]),
]
for x, t, c in methods:
    box(ax, x, 3.7, t, w=2.3, h=1.0, color=c, fc="white", fontsize=9.5, lw=2.4)
    arrow(ax, x, 5.1, x, 4.22, shrink=2, color=c, lw=1.8)

fallback = [
    (1.4, "stuck? try the\ncontrapositive\n$\\neg Q \\Rightarrow \\neg P$"),
    (4.0, "stuck? strengthen\nthe claim, or use\nstrong induction"),
    (6.6, "can't build it?\ncontradiction gives\nexistence, no recipe"),
    (9.2, "derive $R$ and $\\neg R$:\nthe assumption\nwas impossible"),
    (11.7, "none found?\nswitch to proving\nthe $\\forall$ claim"),
]
for x, t in fallback:
    ax.text(x, 2.05, t, ha="center", va="center", fontsize=8.6, color=MUTED)
    arrow(ax, x, 3.18, x, 2.7, shrink=1, lw=1.0)

ax.annotate("", xy=(4.0, 1.35), xytext=(11.7, 1.35),
            arrowprops=dict(arrowstyle="-|>", color=PALETTE[4], lw=1.4, connectionstyle="arc3,rad=0.0", ls="--"))
ax.text(7.9, 0.85, "a failed counterexample hunt often reveals the pattern an induction needs",
        ha="center", fontsize=8.8, color=PALETTE[4])
save(fig, "ch63-proof-strategy")

# ---------------------------------------------------------------- 2. domino picture of induction
fig, ax = diagram(figsize=(10.0, 4.0))
ax.set_xlim(0, 13); ax.set_ylim(-0.6, 4.4)
ax.plot([0.2, 12.8], [0.5, 0.5], color=INK, lw=1.4)


def domino(ax, x, angle_deg, color, label):
    w, h = 0.32, 2.2
    a = np.radians(angle_deg)
    # pivot at the bottom-right corner, tipping to the right
    px, py = x + w, 0.5
    pts = np.array([[-w, 0], [0, 0], [0, h], [-w, h]])
    R = np.array([[np.cos(-a), -np.sin(-a)], [np.sin(-a), np.cos(-a)]])
    p = pts @ R.T + [px, py]
    ax.add_patch(Polygon(p, closed=True, fc=color, ec=INK, lw=1.3, zorder=3))
    ax.text(x + w / 2, 0.12, label, ha="center", va="center", fontsize=10, color=INK)


xs = [0.9 + 1.25 * k for k in range(9)]
angles = [62, 50, 36, 20, 0, 0, 0, 0, 0]
labels = ["$P(0)$", "$P(1)$", "$P(2)$", "$P(3)$", "$P(4)$", "$P(5)$", "$P(6)$", r"$\cdots$", "$P(n)$"]
for k, (x, ang, lab) in enumerate(zip(xs, angles, labels)):
    if lab == r"$\cdots$":
        ax.text(x + 0.16, 1.6, r"$\cdots$", ha="center", va="center", fontsize=18, color=MUTED)
        ax.text(x + 0.16, 0.12, "", ha="center")
        continue
    col = "#fde68a" if ang > 0 else "#dbeafe"
    domino(ax, x, ang, col, lab)

ax.annotate("base case:\npush the first one", xy=(1.25, 1.75), xytext=(0.3, 3.9), fontsize=9.5, color=PALETTE[1],
            ha="left", va="top", arrowprops=dict(arrowstyle="-|>", color=PALETTE[1], lw=1.3))
ax.annotate("", xy=(5.95, 2.5), xytext=(4.95, 2.5),
            arrowprops=dict(arrowstyle="-|>", color=PALETTE[2], lw=1.8, connectionstyle="arc3,rad=-0.45"))
ax.text(5.45, 3.35, "inductive step:\nif $P(k)$ falls, $P(k+1)$ falls", ha="center", fontsize=9.5, color=PALETTE[2])
ax.text(10.6, 3.6, "conclusion: every domino falls,\nhowever far down the line", ha="center", fontsize=9.5, color=INK)
save(fig, "ch63-induction-dominoes")

# ---------------------------------------------------------------- 3. counterexample hunt: circle regions
n = np.arange(1, 11)
regions = np.array([1 + comb(k, 2) + comb(k, 4) for k in n])
doubling = 2.0 ** (n - 1)
fig, ax = new_figure(figsize=(8.0, 4.2))
ax.bar(n - 0.2, doubling, width=0.4, color="#d1d5db", label=r"the pattern you guessed: $2^{n-1}$")
ax.bar(n + 0.2, regions, width=0.4, color=PALETTE[0], label=r"actual regions: $1+\binom{n}{2}+\binom{n}{4}$")
for k, r, d in zip(n, regions, doubling):
    ax.text(k + 0.2, r + 6, str(r), ha="center", fontsize=8.5, color=PALETTE[0])
ax.axvspan(0.5, 5.5, color="#dcfce7", alpha=0.6, zorder=0)
ax.text(3.0, 70, "agree for n = 1..5", ha="center", fontsize=10, color=PALETTE[2])
ax.annotate("first counterexample:\n31, not 32", xy=(6.25, 45), xytext=(5.6, 230), fontsize=10, color=PALETTE[1],
            arrowprops=dict(arrowstyle="-|>", color=PALETTE[1], lw=1.4))
ax.set_xticks(n)
ax.set_xlabel("n = number of points on the circle")
ax.set_ylabel("regions of the disc")
ax.set_ylim(0, 540)
ax.legend(loc="upper left")
ax.set_title("Chords between n points on a circle: a pattern that lies")
save(fig, "ch63-counterexample-hunt")
