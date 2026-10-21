"""Figures for Chapter 50 — Operations Research: Mathematics That Runs the World."""
import sys, pathlib; sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
from style import plt, np, PALETTE, INK, MUTED, new_figure, diagram, node, box, arrow, save

BLUE, RED, GREEN, PURPLE, ORANGE, TEAL, GREY = PALETTE

# ---------------------------------------------------------------- figure 1: two-machine Gantt charts
a = {"A": 4, "B": 7, "C": 2, "D": 6, "E": 3}   # printing times (h)
b = {"A": 5, "B": 2, "C": 6, "D": 3, "E": 7}   # binding times (h)
colors = {"A": BLUE, "B": RED, "C": GREEN, "D": PURPLE, "E": ORANGE}


def schedule(seq):
    t1 = t2 = 0
    out = []
    for j in seq:
        s1 = t1; t1 += a[j]
        s2 = max(t2, t1); t2 = s2 + b[j]
        out.append((j, s1, a[j], s2, b[j]))
    return out, t2


fig, axs = new_figure(2, 1, figsize=(8.6, 4.6), sharex=True)
for ax, seq, title in [(axs[0], "ABCDE", "alphabetical order A-B-C-D-E"),
                       (axs[1], "CEADB", "Johnson's rule C-E-A-D-B")]:
    rows, cmax = schedule(seq)
    for j, s1, d1, s2, d2 in rows:
        for y, s, d in [(1, s1, d1), (0, s2, d2)]:
            ax.broken_barh([(s, d)], (y - 0.35, 0.7), facecolors=colors[j], alpha=0.85,
                           edgecolor="white", lw=1.5)
            ax.text(s + d / 2, y, j, ha="center", va="center", color="white", fontsize=10, weight="bold")
    # idle gaps on the binder
    t = 0
    for j, s1, d1, s2, d2 in rows:
        if s2 > t:
            ax.broken_barh([(t, s2 - t)], (-0.35, 0.7), facecolors="none", edgecolor=MUTED,
                           hatch="///", lw=0.8)
        t = s2 + d2
    ax.axvline(cmax, color=INK, lw=1.4, ls="--")
    ax.text(cmax + 0.3, 0.5, f"makespan\n{cmax} h", va="center", fontsize=9.5, color=INK)
    ax.set_yticks([0, 1]); ax.set_yticklabels(["binder", "printer"])
    ax.set_title(title, fontsize=10.5, loc="left")
    ax.set_ylim(-0.6, 1.6); ax.grid(axis="y", visible=False)
axs[1].set_xlabel("time (hours)")
axs[1].set_xlim(0, 34)
axs[0].text(22.8, 1.0, "hatched = binder idle", fontsize=9, color=MUTED, va="center")
save(fig, "ch50-gantt-johnson")

# ---------------------------------------------------------------- figure 2: EOQ sawtooth and cost curve
D, K, h = 10000.0, 200.0, 4.0
fig, axs = new_figure(1, 2, figsize=(10.2, 3.8))
ax = axs[0]
T = np.linspace(0, 0.3, 3001)
for Q, col, lab in [(1000, BLUE, "Q = 1000 (EOQ)"), (2000, ORANGE, "Q = 2000")]:
    cyc = Q / D
    inv = Q - D * np.mod(T, cyc)
    ax.plot(T * 365, inv, color=col, lw=1.8, label=lab)
    ax.axhline(Q / 2, color=col, lw=1, ls=":")
ax.text(103, 520, "average Q/2", color=MUTED, fontsize=8.5, ha="right")
ax.set_xlabel("day of the year"); ax.set_ylabel("units on the shelf")
ax.set_title("inventory: the sawtooth", fontsize=11)
ax.set_xlim(0, 109.5); ax.set_ylim(0, 2300)
ax.legend(loc="upper right", fontsize=9)

ax = axs[1]
Q = np.linspace(150, 3000, 600)
hold = h * Q / 2; order = K * D / Q
ax.plot(Q, hold, color=GREEN, lw=1.6, label="holding  hQ/2")
ax.plot(Q, order, color=RED, lw=1.6, label="ordering  KD/Q")
ax.plot(Q, hold + order, color=INK, lw=2.4, label="total")
ax.plot(1000, 4000, "o", color=BLUE, ms=7, zorder=5)
ax.annotate("Q* = 1000\ncost 4000 / yr", xy=(1000, 4000), xytext=(1500, 8200), fontsize=9,
            arrowprops=dict(arrowstyle="-|>", color=MUTED), color=INK)
ax.axvspan(667, 1500, color=BLUE, alpha=0.08)
ax.text(1083, 600, "within 8.3%", color=BLUE, fontsize=8.5, ha="center")
ax.set_xlabel("order quantity Q (units)"); ax.set_ylabel("cost per year ($)")
ax.set_title("cost: a flat-bottomed valley", fontsize=11)
ax.set_xlim(0, 3000); ax.set_ylim(0, 12000)
ax.legend(loc="upper right", fontsize=9)
save(fig, "ch50-eoq")
