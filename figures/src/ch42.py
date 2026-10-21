import sys, pathlib; sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
from style import plt, np, PALETTE, INK, MUTED, new_figure, diagram, node, box, arrow, save

# ---------------------------------------------------------------- mechanism as a box
fig, ax = diagram(figsize=(8.0, 4.4))
ax.set_xlim(0, 11)
ax.set_ylim(0, 6)

# agents with private types
for i, (y, lab) in enumerate([(4.8, r"agent 1" "\n" r"type $\theta_1$"),
                              (3.0, r"agent 2" "\n" r"type $\theta_2$"),
                              (1.2, r"agent $n$" "\n" r"type $\theta_n$")]):
    box(ax, 1.2, y, lab, w=1.9, h=1.0, color=PALETTE[i], fontsize=9)
ax.text(1.2, 2.1, r"$\vdots$", ha="center", va="center", fontsize=14, color=MUTED)

# mechanism box
from matplotlib.patches import FancyBboxPatch
ax.add_patch(FancyBboxPatch((4.0, 1.0), 3.0, 4.0, boxstyle="round,pad=0.02,rounding_size=0.2",
                            fc="#f3f4f6", ec=INK, lw=2.2, zorder=1))
ax.text(5.5, 4.55, "MECHANISM", ha="center", va="center", fontsize=11, fontweight="bold", color=INK)
ax.text(5.5, 3.65, r"allocation rule" "\n" r"$x(\hat\theta_1,\dots,\hat\theta_n)$", ha="center", va="center",
        fontsize=9.5, color=INK)
ax.text(5.5, 2.3, r"payment rule" "\n" r"$p_i(\hat\theta_1,\dots,\hat\theta_n)$", ha="center", va="center",
        fontsize=9.5, color=INK)
ax.plot([4.3, 6.7], [3.0, 3.0], color=MUTED, lw=0.8, ls="--")

# reports
for i, y in enumerate([4.8, 3.0, 1.2]):
    ax.annotate("", xy=(4.0, 3.0 + (y - 3.0) * 0.55), xytext=(2.15, y),
                arrowprops=dict(arrowstyle="-|>", color=PALETTE[i], lw=1.6))
ax.text(3.05, 5.15, r"reports $\hat\theta_i$" "\n" "(maybe lies)", ha="center", va="bottom", fontsize=8.5, color=MUTED)

# outputs
box(ax, 9.3, 4.1, "outcome\n(who gets what,\nwhat gets built)", w=2.4, h=1.25, color=PALETTE[2], fontsize=9)
box(ax, 9.3, 1.9, "payments\n(who pays how much)", w=2.4, h=1.0, color=PALETTE[4], fontsize=9)
ax.annotate("", xy=(8.1, 4.1), xytext=(7.0, 3.6), arrowprops=dict(arrowstyle="-|>", color=MUTED, lw=1.6))
ax.annotate("", xy=(8.1, 1.9), xytext=(7.0, 2.3), arrowprops=dict(arrowstyle="-|>", color=MUTED, lw=1.6))

ax.text(5.5, 0.45, r"agent $i$ gets  $u_i = v_i(x;\theta_i) - p_i$   — the designer chooses $x$ and $p$ so that $\hat\theta_i=\theta_i$ is optimal",
        ha="center", va="center", fontsize=9, color=INK)
save(fig, "ch42-mechanism-box")

# ---------------------------------------------------------------- Condorcet cycle
fig, (axl, ax) = plt.subplots(1, 2, figsize=(9.0, 4.0), constrained_layout=True,
                              gridspec_kw={"width_ratios": [1.0, 1.25]})
axl.axis("off")
axl.set_xlim(0, 3.4)
axl.set_ylim(0, 4)
axl.text(1.7, 3.7, "Three voters' rankings", ha="center", fontsize=11, fontweight="bold", color=INK)
rank = [["X", "Y", "Z"], ["Y", "Z", "X"], ["Z", "X", "Y"]]
cols = {"X": PALETTE[0], "Y": PALETTE[1], "Z": PALETTE[2]}
for j, r in enumerate(rank):
    axl.text(0.6 + j * 1.1, 3.15, f"voter {j + 1}", ha="center", fontsize=9.5, color=MUTED)
    for k, c in enumerate(r):
        y = 2.5 - k * 0.8
        axl.add_patch(plt.Rectangle((0.3 + j * 1.1, y - 0.3), 0.6, 0.6, fc="white", ec=cols[c], lw=2))
        axl.text(0.6 + j * 1.1, y, c, ha="center", va="center", fontsize=12, color=cols[c], fontweight="bold")
for k, lab in enumerate(["1st", "2nd", "3rd"]):
    axl.text(0.05, 2.5 - k * 0.8, lab, ha="left", va="center", fontsize=8, color=MUTED)

ax.axis("off")
ax.set_xlim(-1.6, 1.6)
ax.set_ylim(-1.35, 1.5)
ax.set_aspect("equal")
ang = {"X": 90, "Y": 210, "Z": 330}
pos = {k: (np.cos(np.radians(a)), np.sin(np.radians(a))) for k, a in ang.items()}
for k, (x, y) in pos.items():
    ax.add_patch(plt.Circle((x, y), 0.26, fc="white", ec=cols[k], lw=2.4, zorder=3))
    ax.text(x, y, k, ha="center", va="center", fontsize=14, fontweight="bold", color=cols[k], zorder=4)
for a, b in [("X", "Y"), ("Y", "Z"), ("Z", "X")]:
    (x0, y0), (x1, y1) = pos[a], pos[b]
    ax.annotate("", xy=(x1, y1), xytext=(x0, y0),
                arrowprops=dict(arrowstyle="-|>", color=INK, lw=1.8, shrinkA=22, shrinkB=22,
                                connectionstyle="arc3,rad=0.25"))
    mx, my = (x0 + x1) / 2, (y0 + y1) / 2
    r = np.hypot(mx, my)
    ax.text(mx * 1.75 / r * 0.75, my * 1.75 / r * 0.75, f"{a} beats {b}\n2 : 1", ha="center", va="center",
            fontsize=9, color=MUTED)
ax.text(0, 1.45, "Majority preference goes round in a circle", ha="center", fontsize=11, fontweight="bold", color=INK)
save(fig, "ch42-condorcet-cycle")
