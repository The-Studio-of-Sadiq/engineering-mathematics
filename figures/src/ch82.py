import sys, pathlib; sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
from style import plt, np, PALETTE, INK, MUTED, new_figure, diagram, node, box, arrow, save
from matplotlib.patches import Wedge, Rectangle, Circle, FancyArrowPatch

# --- Figure 1: the Ali Baba cave ---------------------------------------------------------
fig, ax = diagram(figsize=(9.6, 5.0))
ax.set_xlim(0, 12); ax.set_ylim(0, 6.6)
cx, cy, R, w = 4.2, 3.55, 2.3, 0.75
ax.add_patch(Wedge((cx, cy), R + w / 2, 0, 360, width=w, fc="#e7dcc8", ec="#8b7355", lw=1.6))
# the entrance tunnel (drawn over the ring's outer wall so the cave is open)
yo = cy - np.sqrt((R + w / 2) ** 2 - (w / 2) ** 2)
ax.add_patch(Rectangle((cx - w / 2 + 0.01, 0.15), w - 0.02, yo + 0.12 - 0.15, fc="#e7dcc8", ec="none", zorder=2))
for sx in (-1, 1):
    ax.plot([cx + sx * w / 2] * 2, [0.15, yo], color="#8b7355", lw=1.6, zorder=3)
# the magic door at the top
ax.add_patch(Rectangle((cx - 0.09, cy + R - w / 2 - 0.02), 0.18, w + 0.04, fc=PALETTE[1], ec=INK, lw=1.2, zorder=4))
ax.text(cx, cy + R + w / 2 + 0.22, "magic door (opens with the secret word)", ha="center", fontsize=9.5, color=PALETTE[1])
# path labels
ax.text(cx - R - 0.85, cy, "path A", ha="center", va="center", fontsize=11, weight="bold", color=PALETTE[0])
ax.text(cx + R + 0.85, cy, "path B", ha="center", va="center", fontsize=11, weight="bold", color=PALETTE[2])
# Peggy deep in path A, Victor at the entrance
pa = (cx + R * np.cos(np.radians(150)), cy + R * np.sin(np.radians(150)))
node(ax, *pa, "P", color=PALETTE[3], r=0.3, fontsize=10)
node(ax, cx, 0.75, "V", color=PALETTE[4], r=0.3, fontsize=10)
ax.text(cx + 0.55, 0.75, "Victor waits at the entrance", va="center", fontsize=9.2, color=INK)
ax.text(pa[0] - 0.15, pa[1] + 0.5, "Peggy picks a path\nin secret", ha="center", fontsize=9.0, color=PALETTE[3])
ax.text(cx, cy, "solid\nrock", ha="center", va="center", fontsize=9, color=MUTED)

# the three steps on the right
steps = [("1  commit", "Peggy walks in by path A or B;\nVictor cannot see which."),
         ("2  challenge", "Victor shouts a random path:\n\"come out by A\" or \"by B\"."),
         ("3  response", "Peggy appears from that path,\nusing the door if she must.")]
for i, (h, t) in enumerate(steps):
    y = 5.3 - 1.55 * i
    ax.text(8.1, y, h, fontsize=10.5, weight="bold", color=INK, va="center")
    ax.text(8.1, y - 0.55, t, fontsize=9.2, color=INK, va="center")
ax.text(8.1, 0.55, "Without the word she survives a round\nwith probability 1/2; after 30 rounds, $10^{-9}$.",
        fontsize=9.0, color=MUTED, va="center")
save(fig, "ch82-ali-baba-cave")

# --- Figure 2: the Schnorr protocol as a message sequence --------------------------------
fig, ax = diagram(figsize=(9.6, 5.4))
ax.set_xlim(0, 12); ax.set_ylim(0, 7.0)
xp, xv = 2.6, 9.4
box(ax, xp, 6.4, "Prover (Peggy)\nknows $x$ with $y=g^x$", w=3.2, h=0.9, color=PALETTE[3], fontsize=10)
box(ax, xv, 6.4, "Verifier (Victor)\nknows $g$, $y$", w=3.2, h=0.9, color=PALETTE[4], fontsize=10)
for xx in (xp, xv):
    ax.plot([xx, xx], [1.5, 5.9], color=MUTED, lw=1.2, ls="--", zorder=1)


def msg(y, left_to_right, label, sub):
    x0, x1 = (xp, xv) if left_to_right else (xv, xp)
    ax.add_patch(FancyArrowPatch((x0, y), (x1, y), arrowstyle="-|>", mutation_scale=16, color=INK, lw=1.6))
    ax.text(6.0, y + 0.18, label, ha="center", va="bottom", fontsize=11, color=INK)
    ax.text(6.0, y - 0.2, sub, ha="center", va="top", fontsize=8.8, color=MUTED)


ax.text(xp - 0.2, 5.15, "pick random $r$\n$t=g^r$", ha="right", va="center", fontsize=9.5, color=INK)
msg(4.85, True, "commitment $t$", "commit")
ax.text(xv + 0.2, 3.85, "pick random\nchallenge $c$", ha="left", va="center", fontsize=9.5, color=INK)
msg(3.45, False, "challenge $c$", "challenge")
ax.text(xp - 0.2, 2.45, "$s=r+cx$\n(mod $q$)", ha="right", va="center", fontsize=9.5, color=INK)
msg(2.05, True, "response $s$", "response")
box(ax, xv, 1.0, "accept iff\n$g^s = t\\,y^c$", w=2.4, h=0.9, color=PALETTE[2], fontsize=10)
ax.text(6.0, 0.3, "$g^s=g^{r+cx}=g^r\\,(g^x)^c=t\\,y^c$: an honest prover always passes",
        ha="center", fontsize=9.2, color=MUTED)
save(fig, "ch82-schnorr-protocol")
