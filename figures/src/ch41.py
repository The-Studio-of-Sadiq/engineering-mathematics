"""Chapter 41 figures: a job-market signaling game tree with information sets; the lemons market unravelling."""
import sys, pathlib; sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
from style import plt, np, PALETTE, INK, MUTED, new_figure, diagram, node, box, arrow, save

BLUE, RED, GREEN, PURPLE, ORANGE, GREY = PALETTE[0], PALETTE[1], PALETTE[2], PALETTE[3], PALETTE[4], PALETTE[6]

# ---------- Figure 1: signaling tree ----------
fig, ax = plt.subplots(figsize=(9.6, 5.6), constrained_layout=True)
ax.set_xlim(0, 12); ax.set_ylim(-0.5, 7); ax.set_aspect("equal"); ax.axis("off")


def seg(x0, y0, x1, y1, label=None, color=INK, off=(0, 0.15), fs=9.5):
    ax.plot([x0, x1], [y0, y1], color=color, lw=1.6, zorder=1)
    if label:
        ax.text((x0 + x1) / 2 + off[0], (y0 + y1) / 2 + off[1], label, ha="center", va="bottom", fontsize=fs,
                color=color)


N = (6, 3.5)
H, L = (6, 5.6), (6, 1.4)
node(ax, *N, "Nature", color=GREY, r=0.55, fontsize=8.5)
seg(*N, *H, "High  (1/3)", color=GREY, off=(0.75, -0.2))
seg(*N, *L, "Low  (2/3)", color=GREY, off=(0.7, -0.2))
node(ax, *H, "W$_H$", color=BLUE, r=0.38)
node(ax, *L, "W$_L$", color=BLUE, r=0.38)
emp = {}
for (wx, wy), t in [(H, "H"), (L, "L")]:
    for side, act in [(-1, "E"), (1, "N")]:
        ex = wx + side * 2.6
        seg(wx, wy, ex, wy, "educate" if act == "E" else "no degree", color=BLUE)
        emp[(t, act)] = (ex, wy)
# information sets: employer sees only the signal
for act in "EN":
    (x, y1), (_, y2) = emp[("H", act)], emp[("L", act)]
    ax.plot([x, x], [y1 - 0.35, y2 + 0.35], ls=(0, (4, 3)), color=RED, lw=1.8, zorder=1)
    ax.text(x + (0.18 if act == "N" else -0.18), 3.5, f"employer sees {act}", rotation=90, ha="center", va="center",
            fontsize=9, color=RED, bbox=dict(fc="white", ec="none", pad=1.5))
for k, v in emp.items():
    node(ax, *v, "Emp", color=RED, r=0.36, fontsize=8)
pay = {("H", "E", "M"): (7, 1), ("H", "E", "C"): (1, 0), ("H", "N", "M"): (10, 1), ("H", "N", "C"): (4, 0),
       ("L", "E", "M"): (3, 0), ("L", "E", "C"): (-3, 1), ("L", "N", "M"): (10, 0), ("L", "N", "C"): (4, 1)}
for (t, act), (x, y) in emp.items():
    side = -1 if act == "E" else 1
    for dy, job in [(0.75, "M"), (-0.75, "C")]:
        tx, ty = x + side * 1.7, y + dy
        seg(x, y, tx, ty, "manager" if job == "M" else "clerk", color=RED, off=(-side * 0.15, 0.24 if dy > 0 else -0.66), fs=8.5)
        w, e = pay[(t, act, job)]
        ax.text(tx + side * 0.15, ty, f"({w}, {e})", ha="left" if side > 0 else "right", va="center", fontsize=9.5,
                color=INK)
ax.text(6, -0.35, "payoffs: (worker, employer).  Manager job worth 10, clerk 4;  degree costs 3 (High) or 7 (Low)",
        ha="center", fontsize=9, color=MUTED)
save(fig, "ch41-signaling-tree")

# ---------- Figure 2: lemons ----------
fig, ax = new_figure(figsize=(6.6, 4.4))
p = np.linspace(0, 1, 200)
ax.plot(p, p, color=INK, lw=2, label="price $p$ asked")
for k, c in [(1.5, RED), (2.5, GREEN)]:
    ax.plot(p, k * p / 2, color=c, lw=2.2, label=f"buyer's value of an average offered car, $k={k}$")
ax.fill_between(p, 0.75 * p, p, color=RED, alpha=0.08)
ax.axhline(0.75, color=GREY, ls=":", lw=1.4)
ax.text(0.02, 0.77, "value of an average car if quality were visible ($k=1.5$)", fontsize=8.5, color=GREY, va="bottom")
ax.annotate("every price is too high:\nbuyers walk away", xy=(0.7, 0.6), xytext=(0.55, 0.22), fontsize=9.5,
            color=RED, arrowprops=dict(arrowstyle="-|>", color=RED, lw=1.2))
ax.set_xlim(0, 1); ax.set_ylim(0, 1.3)
ax.set_xlabel("price $p$  (only sellers with quality $q\\leq p$ offer their car)")
ax.set_ylabel("value")
ax.legend(loc="upper left", fontsize=8.5)
ax.set_title("The market for lemons")
save(fig, "ch41-lemons")
