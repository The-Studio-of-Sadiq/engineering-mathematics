"""Chapter 37 figures: causal junctions, Simpson's paradox, graph surgery."""
import sys, pathlib; sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
from style import plt, np, PALETTE, INK, MUTED, new_figure, diagram, node, box, arrow, save

BLUE, RED, GREEN, PURPLE, ORANGE = PALETTE[:5]

# ---------------------------------------------------------------- 1. three junctions
fig, axes = plt.subplots(1, 3, figsize=(10.5, 3.9), constrained_layout=True)
titles = ["Fork: Z is a confounder", "Chain: M is a mediator", "Collider: C is a common effect"]
notes = ["X and Y correlated through Z.\nConditioning on Z blocks the path:\nADJUST for it.",
         "X affects Y through M.\nConditioning on M blocks the effect:\nDO NOT adjust (for total effect).",
         "X and Y independent.\nConditioning on C opens the path:\nNEVER adjust for it."]
cols = [RED, GREEN, PURPLE]
for k, ax in enumerate(axes):
    ax.set_xlim(0, 6); ax.set_ylim(0, 5.2); ax.set_aspect("equal"); ax.axis("off")
    ax.set_title(titles[k], fontsize=11.5, color=cols[k])
    xs, ys = (1, 2.2), (5, 2.2)
    mid = (3, 4.2) if k != 2 else (3, 0.9)
    lab = ["Z", "M", "C"][k]
    node(ax, *xs, "X", color=BLUE, r=0.5, fontsize=13)
    node(ax, *ys, "Y", color=BLUE, r=0.5, fontsize=13)
    node(ax, *mid, lab, color=cols[k], r=0.5, fontsize=13)
    if k == 0:
        arrow(ax, *mid, *xs, color=INK, shrink=20); arrow(ax, *mid, *ys, color=INK, shrink=20)
    elif k == 1:
        arrow(ax, *xs, *mid, color=INK, shrink=20); arrow(ax, *mid, *ys, color=INK, shrink=20)
    else:
        arrow(ax, *xs, *mid, color=INK, shrink=20); arrow(ax, *ys, *mid, color=INK, shrink=20)
    ny = 0.35 if k != 2 else 4.6
    ax.text(3, ny, notes[k], ha="center", va="center", fontsize=9, color=MUTED)
save(fig, "ch37-junctions")

# ---------------------------------------------------------------- 2. Simpson scatter
rng = np.random.default_rng(37)
fig, ax = new_figure(figsize=(7.2, 4.4))
allx, ally = [], []
groups = [("age 20-35", 1.5, 8.0), ("age 35-50", 3.5, 10.0), ("age 50-65", 5.5, 12.0)]
for (name, cx, cy), c in zip(groups, [BLUE, GREEN, ORANGE]):
    x = cx + rng.normal(0, 0.55, 40)
    y = cy - 0.8 * (x - cx) + rng.normal(0, 0.45, 40)
    ax.scatter(x, y, s=18, color=c, alpha=0.8, label=name, zorder=3)
    xx = np.linspace(x.min(), x.max(), 2)
    b = np.polyfit(x, y, 1)
    ax.plot(xx, np.polyval(b, xx), color=c, lw=2.2, zorder=4)
    allx.append(x); ally.append(y)
allx, ally = np.concatenate(allx), np.concatenate(ally)
b = np.polyfit(allx, ally, 1)
xx = np.linspace(allx.min(), allx.max(), 2)
ax.plot(xx, np.polyval(b, xx), color=RED, lw=2.4, ls="--", label=f"pooled fit (slope {b[0]:+.2f})", zorder=5)
ax.set_xlabel("weekly exercise (hours)")
ax.set_ylabel("cholesterol (arbitrary units)")
ax.set_title("Simpson's paradox: every group slopes down, the pooled data slope up")
ax.legend(loc="upper left", fontsize=9)
save(fig, "ch37-simpson")

# ---------------------------------------------------------------- 3. graph surgery
fig, ax = diagram(figsize=(9.0, 3.6))
ax.set_xlim(0, 12); ax.set_ylim(0, 4.6)
for ox, title in [(0, "Observing:  P(y | x)"), (6.6, "Intervening:  P(y | do(x))")]:
    Z, X, Y = (ox + 2.7, 3.4), (ox + 0.9, 1.2), (ox + 4.5, 1.2)
    ax.text(ox + 2.7, 4.35, title, ha="center", fontsize=11.5, color=INK, weight="bold")
    node(ax, *Z, "Z", color=RED, r=0.45, fontsize=12)
    node(ax, *Y, "Y", color=BLUE, r=0.45, fontsize=12)
    if ox == 0:
        node(ax, *X, "X", color=BLUE, r=0.45, fontsize=12)
        arrow(ax, *Z, *X, color=INK, shrink=30)
    else:
        node(ax, *X, "x", color=ORANGE, r=0.45, fontsize=12, fc="#fff1e6")
        ax.plot([Z[0], X[0]], [Z[1], X[1]], color=MUTED, lw=1.4, ls=":", zorder=1)
        mx, my = (Z[0] + X[0]) / 2, (Z[1] + X[1]) / 2
        ax.plot([mx - 0.2, mx + 0.2], [my - 0.2, my + 0.2], color=RED, lw=3)
        ax.plot([mx - 0.2, mx + 0.2], [my + 0.2, my - 0.2], color=RED, lw=3)
        ax.text(mx - 0.35, my + 0.1, "arrow cut", ha="right", fontsize=9, color=RED)
        ax.text(X[0], 0.3, "set by hand", ha="center", fontsize=9, color=ORANGE)
    arrow(ax, *Z, *Y, color=INK, shrink=30)
    arrow(ax, *X, *Y, color=INK, shrink=30)
ax.plot([6.0, 6.0], [0.2, 4.2], color="#e5e7eb", lw=1.5)
save(fig, "ch37-surgery")
