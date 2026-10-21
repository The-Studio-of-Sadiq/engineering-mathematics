import sys, pathlib; sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
from style import plt, np, PALETTE, INK, MUTED, new_figure, diagram, node, box, arrow, save

# ---------------------------------------------------------------- bid shading
fig, ax = new_figure(figsize=(7.0, 4.4))
v = np.linspace(0, 1, 200)
ax.plot(v, v, color=INK, lw=2.4, label="second-price: bid your value, $b=v$")
for k, n in enumerate([2, 3, 5, 10]):
    ax.plot(v, (n - 1) / n * v, color=PALETTE[k], lw=2,
            label=fr"first-price, $n={n}$: $b=\frac{{{n - 1}}}{{{n}}}v$")
ax.fill_between(v, v / 2, v, color=PALETTE[0], alpha=0.07)
ax.annotate("shading\n(n = 2)", xy=(0.8, 0.6), xytext=(0.88, 0.32), fontsize=9, color=MUTED,
            ha="center", arrowprops=dict(arrowstyle="-", color=MUTED, lw=1))
ax.set_xlim(0, 1)
ax.set_ylim(0, 1.02)
ax.set_xlabel("private value $v$ (values i.i.d. uniform on [0, 1])")
ax.set_ylabel("equilibrium bid $b(v)$")
ax.set_title("More rivals, less shading")
ax.legend(loc="upper left", fontsize=9)
save(fig, "ch43-bid-shading")

# ---------------------------------------------------------------- package auction allocation
fig, ax = diagram(figsize=(9.0, 4.4))
ax.set_xlim(0, 12)
ax.set_ylim(0, 6)
win = {"A": PALETTE[0], "B": PALETTE[0], "C": PALETTE[2], "D": PALETTE[4]}
lic = {"A": (0.6, 3.0), "B": (3.0, 3.0), "C": (0.6, 0.6), "D": (3.0, 0.6)}
names = {"A": "North-west", "B": "North-east", "C": "South-west", "D": "South-east"}
for k, (x, y) in lic.items():
    ax.add_patch(plt.Rectangle((x, y), 2.3, 2.3, fc=win[k], ec="white", lw=3, alpha=0.22))
    ax.add_patch(plt.Rectangle((x, y), 2.3, 2.3, fc="none", ec=win[k], lw=2))
    ax.text(x + 1.15, y + 1.35, k, ha="center", va="center", fontsize=20, fontweight="bold", color=win[k])
    ax.text(x + 1.15, y + 0.6, names[k] + " licence", ha="center", va="center", fontsize=8.5, color=INK)
ax.text(2.95, 5.65, "Winning allocation", ha="center", fontsize=11, fontweight="bold", color=INK)

rows = [("Bidder 1", "{A, B}", 10, True, PALETTE[0]),
        ("Bidder 2", "{C}", 6, True, PALETTE[2]),
        ("Bidder 3", "{D}", 5, True, PALETTE[4]),
        ("Bidder 4", "{B, C}", 12, False, MUTED),
        ("National", "{A, B, C, D}", 20, False, MUTED)]
ax.text(8.9, 5.65, "Package bids", ha="center", fontsize=11, fontweight="bold", color=INK)
for i, (who, pkg, val, won, c) in enumerate(rows):
    y = 4.9 - i * 0.78
    ax.text(6.4, y, who, ha="left", va="center", fontsize=10, color=c, fontweight="bold" if won else "normal")
    ax.text(8.2, y, pkg, ha="left", va="center", fontsize=10, color=INK)
    ax.text(10.6, y, f"{val}", ha="right", va="center", fontsize=10, color=INK)
    ax.text(11.0, y, "wins" if won else "loses", ha="left", va="center", fontsize=9, color=c)
ax.plot([6.3, 11.8], [1.0, 1.0], color=MUTED, lw=0.8)
ax.text(6.4, 0.6, "best split: 10 + 6 + 5 = 21  >  20 (national)  >  12 + 5 = 17",
        ha="left", va="center", fontsize=9, color=INK)
save(fig, "ch43-spectrum-allocation")
